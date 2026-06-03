library(cellgpsr)
library(readr)
library(dplyr)

cells_csv <- Sys.getenv("CELLGPSR_LN_CELLS_CSV", unset = "cells.csv.gz")
celltypes_csv <- Sys.getenv("CELLGPSR_LN_CELLTYPES_CSV", unset = "cell_types.csv")
transcripts_parquet <- Sys.getenv("CELLGPSR_LN_TRANSCRIPTS_PARQUET", unset = "transcripts.parquet")
output_dir <- Sys.getenv("CELLGPSR_OUTPUT_DIR", unset = ".")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(cells_csv, celltypes_csv, transcripts_parquet)) {
  if (!file.exists(path)) {
    stop("Set CELLGPSR_LN_CELLS_CSV, CELLGPSR_LN_CELLTYPES_CSV and CELLGPSR_LN_TRANSCRIPTS_PARQUET before running this template.")
  }
}

cells = read_csv(cells_csv)

celltype <- read_csv(celltypes_csv)

cells <- cells %>%
  left_join(celltype %>% select(cell_id, group), by = "cell_id")
cells$group <- paste0("$\\textbf{", gsub(" ", "~", cells$group), "}$")

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  cells,
  cluster_col = "group",
  x_col = "x_centroid",
  y_col = "y_centroid",
  method = "average"
)

# 绘制 cophenetic 热图
plot_cophenetic_heatmap(
  result[['row_cophenetic_df']],
  figsize = c(15, 15),
  cellwidth = 13,
  matrix_name = "row_coph",
  output_dir = output_dir,
  sample = "LN"
)


# 读取transcript ------------------------------------------------------------

library(arrow)
library(dplyr)

# 打开 Parquet 数据集（不会立刻加载所有数据）
ds <- open_dataset(transcripts_parquet, format = "parquet")

# 通过 dplyr 操作进行筛选或选择部分列，例如只读取部分列和符合条件的行
freq_table <- ds %>%
  count(feature_name) %>%   # 计算每个 feature_name 的出现次数
  collect()                 # 收集结果到内存

freq_table = subset(freq_table, grepl("^(CCR|CCL|CXCR|CXCL)", feature_name))

subset_ds <- ds %>%
  select(x_location, y_location, z_location, cell_id, feature_name) %>%
  semi_join(freq_table, by = "feature_name")

df <- collect(subset_ds)
df$feature_name <- paste0("$\\textit{", df$feature_name, "}$")


# 合并cell和transcript的StructureMap ------------------------------------------

matrix <- rbind(setNames(cells[, c("x_centroid","y_centroid","cell_id","group")], c("x","y","cell_id","group")),
                setNames(df[, c("x_location","y_location","cell_id","feature_name")], c("x","y","cell_id","group")))

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  matrix,
  cluster_col = "group",
  x_col = "x",
  y_col = "y",
  method = "average"
)

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

# 绘制 cophenetic 热图
plot_cophenetic_heatmap(
  result[['row_cophenetic_df']],
  figsize = c(15, 15),
  matrix_name = "row_coph",
  output_dir = output_dir,
  sample = "5K_LN"
)
