library(cellgpsr)
library(readr)

transcripts_csv <- Sys.getenv("CELLGPSR_NPC_TRANSCRIPTS_CSV", unset = "874141_transcripts.csv")
cells_csv <- Sys.getenv("CELLGPSR_NPC_CELLS_CSV", unset = "874141_celltype.csv")
if (!file.exists(transcripts_csv)) {
  stop("Set CELLGPSR_NPC_TRANSCRIPTS_CSV to the NPC transcript CSV file.")
}
if (!file.exists(cells_csv)) {
  stop("Set CELLGPSR_NPC_CELLS_CSV to the NPC cell-type CSV file.")
}
X874141_transcripts <- read_csv(transcripts_csv)


# 合并cell和transcript的StructureMap ------------------------------------------

# 读取 CSV 文件
transcript = subset(X874141_transcripts, grepl("^CXCL5", feature_name))
transcript$feature_name <- paste0("$\\textit{", transcript$feature_name, "}$")

cells = read.csv(cells_csv)
cells$feature_name <- paste0("$\\textbf{", gsub(" ", "~", cells$feature_name), "}$")

df = rbind(transcript, cells)

# 计算 cophenetic 距离
result <- compute_cophenetic_distances_from_df(
  df,
  cluster_col = "feature_name",
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
  output_dir = './',
  sample = "NPC_CXCL5"
)

