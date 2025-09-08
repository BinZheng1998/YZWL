setwd("~/project/03_3D_spatial/04_cell_annotation/new_data_20250827/")
library(Seurat)
library(dplyr)

####E2.5
data <- readRDS("E25_gut.rds")
head(data@meta.data,n=100)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(2,3)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
head(data@meta.data)
data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Proventriculus Progenitor Cells",
  "8" = "Proventriculus Progenitor Cells",
  "9" = "Proventriculus Progenitor Cells",
  "14" = "Proventriculus Progenitor Cells",
  "1" = "Foregut Progenitor Cells1",
  "2" = "Hepatocyte Progenitor Cells",
  "24" = "Hepatocyte Progenitor Cells",
  "3" = "Duodenum-Midgut Progenitor Cells",
  "20" = "Duodenum-Midgut Progenitor Cells",
  "4" = "Gizzard Progenitor Cells",
  "11" = "Gizzard Progenitor Cells",
  "15" = "Gizzard Progenitor Cells",
  "18" = "Gizzard Progenitor Cells",
  "19" = "Gizzard Progenitor Cells",
  "22" = "Gizzard Progenitor Cells",
  "23" = "Gizzard Progenitor Cells",
  "25" = "Gizzard Progenitor Cells",
  "28" = "Gizzard Progenitor Cells",
  "5" = "Lung Bud1",
  "6" = "Lung Bud1",
  "17" = "Lung Bud1",
  "26" = "Lung Bud2",
  "7" = "Foregut Progenitor Cells2",
  "31" = "Foregut Progenitor Cells2",
  "32" = "Foregut Progenitor Cells2",
  "10" = "Gizzard-Duodenum Progenitor Cells",
  "13" = "Foregut Progenitor Cells3",
  "29" = "Foregut Progenitor Cells3",
  "16" = "Hepatocyte Mesenchymal Progenitor Cells",
  "21" = "Pancreatic Bud",
  "27" = "Duodenum Progenitor Cells",
  "34" = "Duodenum Progenitor Cells",
  "30" = "Hepatopancreatic Ductal Progenitor Cells",
  "33" = "Foregut Mucosa Epithelial Progenitor Cells"
)
head(data@meta.data,n=50)
table(data@meta.data$celltype_gut)

cluster <- data.frame(data@meta.data$Spatial_snn_res.3,data@meta.data$celltype_gut)
cluster1 <- unique(cluster)
cell_number <- data.frame(table(data@meta.data$Spatial_snn_res.3))
colnames(cell_number) <- c("Cluster",'Cell_number')
cell_number$Prop <- cell_number$Cell_number/sum(cell_number$Cell_number)
cell_number$Cell <- cluster1$data.meta.data.celltype_gut[match(cell_number$Cluster,cluster1$data.meta.data.Spatial_snn_res.3)]
write.table(cell_number,'E2.5_cell_number.txt',sep = '\t',row.names = F)

#####E3.5
data <- readRDS("E35_gut.rds")
head(data@meta.data,n=100)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(2,3)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
head(data@meta.data)

data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Duodenum-Midgut Progenitor Cells",
  "1" = "Hindgut Progenitor Cells",
  "8" = "Hindgut Progenitor Cells",
  "29" = "Hindgut Progenitor Cells",
  "3" = "Gizzard Progenitor Cells",
  "11" = "Gizzard Progenitor Cells",
  "16" = "Gizzard Progenitor Cells",
  "21" = "Gizzard Progenitor Cells",
  "26" = "Gizzard Progenitor Cells",
  "4" = "Proventriculus Progenitor Cells",
  "5" = "Proventriculus Progenitor Cells",
  "18" = "Proventriculus Progenitor Cells",
  "20" = "Proventriculus Progenitor Cells",
  "6" = "Duodenum Enterocyte-Neurons Cells",
  "15" = "Duodenum Enterocyte-Neurons Cells",
  "7" = "Hepatocyte Progenitor Cells",
  "10" = "Hepatocyte Progenitor Cells",
  "24" = "Hepatocyte Progenitor Cells",
  "32" = "Hepatocyte Progenitor Cells",
  "34" = "Hepatocyte Progenitor Cells",
  "38" = "Hepatocyte Progenitor Cells",
  "39" = "Hepatocyte Progenitor Cells",
  "9" = "Gizzard-Duodenum Progenitor Cells",
  "12" = "Gizzard-Duodenum Progenitor Cells",
  "13" = "Foregut Mesothelial cells",
  "25" = "Foregut Mesothelial cells",
  "28" = "Foregut Mesothelial cells",
  "36" = "Foregut Mesothelial cells",
  "14" = "Duodenum Progenitor Cells",
  "17" = "Duodenum Progenitor Cells",
  "19" = "Hindgut Mucosa Epithelial Cells",
  "22" = "Lung Mesenchymal Progenitor Cells",
  "27" = "Lung Mesenchymal Progenitor Cells",
  "23" = "Lung Secondary Bronchi",
  "33" = "Lung Primary Bronchi",
  "30" = "Colacal Epithelial Cells",
  "35" = "Proventriculus Mucosa Epithelial Cells",
  "37" = "Gizzard Mucosa Epithelial Cells",
  "40" = "Duodenum Mucosa Epithelial Cells",
  "41" = "Pancreatic Bud",
  "42" = "Pancreatic Bud"
)

head(data@meta.data,n=50)
table(data@meta.data$celltype_gut)

cluster <- data.frame(data@meta.data$Spatial_snn_res.3,data@meta.data$celltype_gut)
cluster1 <- unique(cluster)
cell_number <- data.frame(table(data@meta.data$Spatial_snn_res.3))
colnames(cell_number) <- c("Cluster",'Cell_number')
cell_number$Prop <- cell_number$Cell_number/sum(cell_number$Cell_number)
cell_number$Cell <- cluster1$data.meta.data.celltype_gut[match(cell_number$Cluster,cluster1$data.meta.data.Spatial_snn_res.3)]
write.table(cell_number,'E3.5_cell_number.txt',sep = '\t',row.names = F)


####E4.5
data <- readRDS('E45_gut.rds')
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(1,2,3)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
head(data@meta.data)
data@meta.data$celltype_gut <- recode(
  data@meta.data$celltype_gut,
  "0" = "Hepatocyte Progenitor Cells",
  "2" = "Hepatocyte Progenitor Cells",
  "5" = "Hepatocyte Progenitor Cells",
  "31" = "Hepatocyte Progenitor Cells",
  "40" = "Hepatocyte Progenitor Cells",
  "1" = "Gizzard Progenitor Cells",
  "21" = "Gizzard Progenitor Cells",
  "23" = "Gizzard Progenitor Cells",
  "41" = "Gizzard Progenitor Cells", 
  "3" = "Gizzard-Duodenum Mesenchymal Progenitor Cells",
  "4" = "Trachea Mesenchymal Progenitor Cells", 
  "6" = "Lung Mesenchymal Progenitor Cells",
  "49" = "Lung Mesenchymal Progenitor Cells",
  "7" = "Midgut-Hindgut Progenitor Cells",
  "17" = "Midgut-Hindgut Progenitor Cells",
  "43" = "Midgut-Hindgut Smooth Muscle Progenitors",
  "8" = "Foregut Mesothelial cells",
  "12" = "Foregut Mesothelial cells",
  "25" = "Foregut Mesothelial cells",
  "29" = "Foregut Mesothelial cells",
  "9" = "Proventriculus Progenitor Cells",
  "15" = "Proventriculus Progenitor Cells",
  "30" = "Proventriculus Progenitor Cells",
  "32" = "Proventriculus Progenitor Cells",
  "42" = "Proventriculus Progenitor Cells",
  "10" = "Lung Progenitor Cells",
  "20" = "Lung Secondary Bronchi",
  "33" = "Lung Primary Bronchi",
  "11" = "Gizzard Mesenchymal Progenitor Cells",
  "19" = "Gizzard Mesenchymal Progenitor Cells",
  "44" = "Gizzard Mesenchymal Progenitor Cells",
  "50" = "Gizzard Mesenchymal Progenitor Cells",
  "13" = "Duodenum-Midgut Progenitor Cells",
  "14" = "Colacal Epithelial Cells",
  "16" = "Hindgut Progenitor Cells",
  "27" = "Hindgut Progenitor Cells",
  "36" = "Hindgut Progenitor Cells",
  "18" = "Hindgut Smooth Muscle Progenitor Cells",
  "22" = "Duodenum Mesenchymal Progenitor Cells",
  "35" = "Duodenum Progenitor Cells",
  "45" = "Duodenum Progenitor Cells",
  "46" = "Duodenum Progenitor Cells",
  "47" = "Duodenum Progenitor Cells",
  "24" = "Duodenum Enterocyte-Neurons Cells",
  "38" = "Duodenum Enterocyte-Neurons Cells",
  "26" = "Hindgut Mucosa Epithelial Cells",
  "28" = "Gizzard Mucosa Epithelial Cells",
  "34" = "Proventriculus Mucosa Epithelial Cells",
  "48" = "Duodenum Mucosa Epithelial Cells",
  "37" = "Pancreatic Bud",
  "52" = "Pancreatic Bud",
  "39" = "Spleno-Pancreatic Mesenchyme Progenitor Cells",
  "51" = "Hepatopancreatic Ductal Progenitor Cells"
)
head(data@meta.data,n=50)
table(data@meta.data$Spatial_snn_res.3)

cluster <- data.frame(data@meta.data$Spatial_snn_res.3,data@meta.data$celltype_gut)
cluster1 <- unique(cluster)
cell_number <- data.frame(table(data@meta.data$Spatial_snn_res.3))
colnames(cell_number) <- c("Cluster",'Cell_number')
cell_number$Prop <- cell_number$Cell_number/sum(cell_number$Cell_number)
cell_number$Cell <- cluster1$data.meta.data.celltype_gut[match(cell_number$Cluster,cluster1$data.meta.data.Spatial_snn_res.3)]
write.table(cell_number,'E4.5_cell_number.txt',sep = '\t',row.names = F)
