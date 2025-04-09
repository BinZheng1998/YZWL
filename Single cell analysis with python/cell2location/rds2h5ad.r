library(Seurat)
library(SeuratDisk)

obj <- readRDS("/share/appspace_data/shared_groups/BGI/04.Project/xuchunyan/Chicken/05.SC_analysis/03.Merge/merge_results/E3.5.merged.filter_norm_rs0.8.rds")

## cell2loc 需要counts矩阵，SaveH5Seurat默认保存scale.data和data，剔除scale.data后则保存data和counts
obj <- DietSeurat(
  obj,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)

SaveH5Seurat(obj,filename="E3.5.merged.filter_norm_rs0.8.h5seurat", overwrite = TRUE)
Convert("E3.5.merged.filter_norm_rs0.8.h5seurat", assay = "RNA", dest = "h5ad", overwrite = TRUE)
