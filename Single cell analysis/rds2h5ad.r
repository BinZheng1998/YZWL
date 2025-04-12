obj <- readRDS("E2.5_C24_C27_res0.3.processed.rds")
head(obj@assays$SCT@counts)
DefaultAssay(obj) <- "SCT"
print(head(obj@assays$SCT@counts[, 1:5]))
obj <- DietSeurat(
  obj,
  counts = TRUE,       # 保留 SCT@counts
  data = FALSE,         # 不导出 SCT@data (标准化后的值)
  scale.data = FALSE,   # 不导出 SCT@scale.data
  assays = "SCT"        # 明确指定处理 SCT assay
)
SaveH5Seurat(
  obj,
  filename = "E3.5.merged.filter_norm_rs0.8.h5seurat",
  overwrite = TRUE,
  assay = "SCT",        # 指定导出 SCT assay
  slots = "counts"      # 强制导出 counts 层
)

Convert(
  "E3.5.merged.filter_norm_rs0.8.h5seurat",
  dest = "h5ad",
  overwrite = TRUE,
  assay = "SCT",        # 指定读取 SCT assay
  slot = "counts"       # 明确导出 counts 层
)
