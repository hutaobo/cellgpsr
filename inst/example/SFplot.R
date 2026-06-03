# 加载必要的 R 包
library(dplyr)
library(ggplot2)

seurat_rds <- Sys.getenv(
  "CELLGPSR_FIBROSIS_SEURAT_RDS",
  unset = "GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds"
)
if (!file.exists(seurat_rds)) {
  stop("Set CELLGPSR_FIBROSIS_SEURAT_RDS to the fibrosis Seurat RDS file.")
}
data <- readRDS(seurat_rds)

metadata <- data@meta.data
samplelist <- unique(metadata$sample)

for (sample in samplelist) {
  submeta <- metadata[metadata$sample == sample, ]
  result <- compute_cluster_average_nn_distance_matrix(submeta,
                                                       cluster_col = "final_CT",
                                                       x_col = "x_centroid",
                                                       y_col = "y_centroid")
  plot_cophenetic_heatmap(
    matrix = result,
    matrix_name = "",
    figsize = c(12, 12),
    output_dir = "./output",
    sample = sample
  )
}
