source('~/project/03_3D_spatial/03-script/export_clusters_to_ply.R')
#去除散点
#E25
library(dbscan)
E25 <- readRDS('E25_lung_260126.rds')
E25 <- NormalizeData(E25,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
coord <- read.csv('../E2.5_251107_coord_3d.csv')
coord$X <- sub("-[^-]*$", "", coord$X)
E25$trans_x <- NULL
E25$trans_y <- NULL
E25$trans_z <- NULL
colnames(coord)[1] <- 'cell'
rownames(coord) <- coord$cell
E25 <- RenameCells(object = E25, new.names = E25$cell)
E25 <- AddMetaData(E25,coord)
head(E25)
DimPlot(E25mes,group.by = 'lung_celltype')

export_clusters_to_ply(E25,output_dir = './E25mes',coord_cols = c("trans_x", "trans_y", "trans_z"),cluster_col = 'celltype')

coords_3d <- E25@meta.data[, c("trans_x", "trans_y", "trans_z")]
db_res <- dbscan(coords_3d, eps = 20, minPts = 5)
E25$cluster_3d <- db_res$cluster
DimPlot(E25,group.by = 'cluster_3d')
export_clusters_to_ply(E25,output_dir = './E25',coord_cols = c("trans_x", "trans_y", "trans_z"),cluster_col = 'cluster_3d')

#保存散点信息
E25_unknown <- subset(E25, subset = cluster_3d == "0",)
E25_cell <- data.frame(E25_unknown$cell,E25_unknown$cluster_3d)
E25_cell$E25_unknown.cluster_3d[E25_cell$E25_unknown.cluster_3d == '0'] <- 'Unknown'
colnames(E25_cell) <- c('cell','celltype')
write.csv(E25_cell,'E25_lung_unknown_cell_260127.csv')

E25_new <- subset(E25, subset = cluster_3d == "1",)
saveRDS(E25_new,'E25_lung_260127.rds')

#E35
E35 <- readRDS('E35_lung_260126.rds')
E35 <- NormalizeData(E35,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(E35,group.by = 'celltype')

coords_3d <- E35@meta.data[, c("trans_x", "trans_y", "trans_z")]
db_res <- dbscan(coords_3d, eps = 20, minPts = 5)
E35$cluster_3d <- db_res$cluster
DimPlot(E35,group.by = 'cluster_3d')
export_clusters_to_ply(E35,output_dir = './E35',coord_cols = c("trans_x", "trans_y", "trans_z"),cluster_col = 'cluster_3d')

#保存散点信息
E35_unknown <- subset(E35, subset = cluster_3d == "0",)
E35_cell <- data.frame(E35_unknown$cell,E35_unknown$cluster_3d)
E35_cell$E35_unknown.cluster_3d[E35_cell$E35_unknown.cluster_3d == '0'] <- 'Unknown'
colnames(E35_cell) <- c('cell','celltype')
write.csv(E35_cell,'E35_lung_unknown_cell_260127.csv')

E35_new <- subset(E35, subset = cluster_3d != "0",)
DimPlot(E35_new,group.by = 'cluster_3d')
saveRDS(E35_new,'E35_lung_260127.rds')

source('~/project/03_3D_spatial/03-script/export_clusters_to_ply.R')
#E45
E45 <- readRDS('E45_lung_260126.rds')
E45 <- NormalizeData(E45,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(E45,group.by = 'celltype')

coords_3d <- E45@meta.data[, c("trans_x", "trans_y", "trans_z")]
db_res <- dbscan(coords_3d, eps = 20, minPts = 5)
E45$cluster_3d <- db_res$cluster
DimPlot(E45,group.by = 'cluster_3d')
export_clusters_to_ply(E45,output_dir = './E45',coord_cols = c("trans_x", "trans_y", "trans_z"),cluster_col = 'cluster_3d')

#保存散点信息
E45_unknown <- subset(E45, subset = cluster_3d %in% c("0","13",'14'),)
E45_cell <- data.frame(E45_unknown$cell,E45_unknown$cluster_3d)
E45_cell$E45_unknown.cluster_3d[E45_cell$E45_unknown.cluster_3d == '0'] <- 'Unknown'
colnames(E45_cell) <- c('cell','celltype')
write.csv(E45_cell,'E45_lung_unknown_cell_260127.csv')

E45_new <- subset(E45, subset = !cluster_3d %in% c("0","13",'14'),)
E45_new
DimPlot(E45_new,group.by = 'cluster_3d')
head(E45_new)
saveRDS(E45_new,'E45_lung_260127.rds')
