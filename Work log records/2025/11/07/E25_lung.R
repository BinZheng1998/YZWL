setwd('~/project/03_3D_spatial/02_result/251107_E25_lung/')
library(Seurat)
library(scplotter)
df <- readRDS("~/project/03_3D_spatial/00_data/E2.5_SP_new_coord_3d_2511.rds")
df <- UpdateSeuratObject(df)
DefaultAssay(df) <- 'Spatial'

color_mapping <- export_clusters_to_ply(df,
                                        output_dir = "../lung/",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "mapped_celltype")
head(df)
DimPlot(df,group.by = 'mapped_celltype',raster = F,label = T)+NoLegend()
FeatureStatPlot(df,plot_type = 'dim',features = 'NKX2-1')
FeaturePlot(df,features = 'NKX2-1',raster = F)

data <- readRDS('~/project/03_3D_spatial/00_data/lung/E25_new_lung_newcoord.rds')
table(data$celltype_gut)
head(data1)
data1 <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.3,0.5,1,1.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(data1,group.by = 'Spatial_snn_res.1')
FeaturePlot(data1,features = 'AGR2',raster = F)

cellID <- data$cell

df$lung <- 'unknown'
df$lung[df$cell %in% cellID] <- 'lung'
table(df$lung)

DimPlot(df,group.by = 'lung',raster = F)

table(df$mapped_celltype)
df1 <- subset(df,subset = mapped_celltype %in% c('Pharyngeal Epithelium','Pharyngeal Ectoderm'))

DimPlot(df1,group.by = 'mapped_celltype',raster = F,reduction = 'umap')
FeaturePlot(df1,features = 'NKX2-1',raster = F)

df2 <- NormalizeData(df1,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5,1,1.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
FeaturePlot(df2,features = 'NKX2-1',raster = F)
DimPlot(df2,group.by = 'Spatial_snn_res.1.5',raster = F,reduction = 'umap',label = T)
DimPlot(df2,group.by = 'lung',raster = F)


df3 <- subset(df2,subset = Spatial_snn_res.1.5 %in% c('1','3','11','18','10','14','16') )
df3 <- NormalizeData(df3,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5,1,1.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(df3,group.by = 'Spatial_snn_res.1',raster = F,reduction = 'umap',label = T)
DimPlot(df3,group.by = 'lung',raster = F)
color_mapping <- export_clusters_to_ply(df3,
                                        output_dir = "../lung/",
                                        coord_cols = c("trans_x", "trans_y", "trans_z"),
                                        cluster_col = "Spatial_snn_res.1.5")
FeaturePlot(df3,features = 'RIPPLY3',raster = F)
all_markers <- FindAllMarkers(object = df3,group.by = 'Spatial_snn_res.1.5',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)

FeaturePlot(df3,features = 'SHH',raster = F)
