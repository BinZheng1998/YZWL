setwd('~/project/03_3D_spatial/02_result/251029_lung_recluster/E4.5/')
library(Seurat)
library(dplyr)

df <- readRDS("../../04_cell_annotation/new_data_20250827/E4.5_foregut.rds")
table(df@meta.data$celltype_gut)
df <- UpdateSeuratObject(df)
df1 <- subset(df,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells',
                                              'Lung Progenitor Cells','Lung Primary Bronchi',
                                              'Lung Secondary Bronchi','Trachea Mesenchymal Progenitor Cells'))

#df1 <- subset(df,subset = celltype_gut %in% c('Lung Primary Bronchi'))
head(df1)

FeaturePlot(data,features = c("FGF10","BMP4"))
DimPlot(df1, reduction = "umap", group.by = c("celltype_gut"),label = T)+NoLegend()

data <- NormalizeData(df1,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.3,0.4,0.5,0.7,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(data, reduction = "umap", group.by = c("Spatial_snn_res.0.3","celltype_gut"),label = T)+NoLegend()
FeaturePlot(data,features = c("FGF10"))

#为每个cluster生成单独的cellID文件
cell_ids <- rownames(data@meta.data)
cell_types <- data@meta.data$Spatial_snn_res.0.3
unique_celltypes <- unique(cell_types)
if (!dir.exists("celltype_ids")) {
  dir.create("celltype_ids")
}
for (ct in unique_celltypes) {
  ct_cell_ids <- cell_ids[cell_types == ct]
  safe_ct_name <- gsub("[^[:alnum:]_]", "_", ct)
  filename <- paste0("celltype_ids/", safe_ct_name, "_cellIDs.txt")
  
  write.table(ct_cell_ids, file = filename, 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat(paste0("Saved ", length(ct_cell_ids), " cells for: ", ct, "\n"))
}

all_markers <- FindAllMarkers(object = data,group.by = 'Spatial_snn_res.0.5',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
