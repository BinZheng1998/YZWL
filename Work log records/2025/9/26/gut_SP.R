setwd(dir = "~/project/03_3D_spatial/04_cell_annotation/new_data_2509/")
library(Seurat)
library(dplyr)
library(ggplot2)
library(readxl)
########################################################################2.5
df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E2.5_new_sp.rds")
df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
df@meta.data
#data1 <- subset(df,subset = Tissue %in% c("Heart"))
df@meta.data$Celltype <- gsub("Pharyngeal Neural Crest-derived Sensory/Motor Progenitors","Pharyngeal Neural Crest-derived Sensory Progenitors",df@meta.data$Celltype)
data <- subset(df,subset = Celltype %in% c("Digestive Gland Progenitors",
                                           "Stomach Mesenchyme","Foregut Primordium",
                                           "Outflow Tract","Hepatocytes","Hepatic Mesothelium",
                                           "Pancreatic Primordium","Periportal Mesenchyme",
                                           "Midgut Primordium","Splanchnic Mesoderm",
                                           "Neural Crest Cells","Migratory Neural Crest Cells"))
data@meta.data
table(df@meta.data$Celltype)
DimPlot(data, reduction = "umap",group.by = "Celltype")

DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,4,4.5,5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
p1<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.1",label = T) + ggtitle("Res 0.1")
p2<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.0.2",label = T) + ggtitle("Res 0.2")
p3<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.2",label = T) + ggtitle("Res 2")
p3
p4<-DimPlot(data, reduction = "umap", group.by = "Spatial_snn_res.3",label = T) + ggtitle("Res 3")
p4
DimPlot(data, reduction = "umap", group.by = "Celltype",label = T) 


########################################################################3.5
df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E3.5_new_sp.rds")
df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
df@meta.data
data <- subset(df,subset = Celltype %in% c("Midgut Mesenchyme","Midgut Epithelium",
                                           "Hindgut Epithelium","Hindgut Mesenchyme",
                                           "Foregut Epithelium",
                                            #"Peritoneal Mesothelium",
                                           "Hepatocytes","Hepatic Mesothelium","Hepatic Primordium",
                                           "Cloacal Epithelium",
                                           "Gizzard Mesenchyme","Proventriculus Mesenchyme",
                                           "Lung Epithelium","Lung Mesenchyme",
                                           "Pancreatic Primordium","Neural Crest Cells"
                                            ))


head(data@meta.data)
DimPlot(data, reduction = "umap", group.by = "Celltype",raster = F,)
DimPlot(df, reduction = "umap", group.by = "Celltype",raster = F,)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(2,3,4,5)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)

########################################################################4.5
df <- readRDS("~/project//03_3D_spatial//04_cell_annotation//new_data_20250827//E4.5_new_sp.rds")
#df@meta.data$Celltype[is.na(df@meta.data$Celltype)] <- "Unidentified"
df@meta.data
data <- subset(df,subset = Celltype %in% c("Gut Smooth Muscle Progenitors","Gastric Epithelium",
                                           "Midgut Smooth Muscle Progenitors","Midgut Epithelium","Midgut Mesenchyme",
                                           #"Gonadal Ridge Mesenchyme","Gonadal Primordium",
                                           "Proventriculus Mesenchyme","Gizzard Mesenchyme",
                                           "Cloacal Epithelium",
                                           "Hepatic Mesenchyme","Hepatic Endoderm","Hepatocytes",
                                           "Foregut Mesenchyme","Foregut Epithelium",
                                           "Hindgut Smooth Muscle Progenitors","Hindgut Epithelium",
                                           "Pancreatic Primordium",
                                           "Lung Epithelium","Lung Mesenchyme",
                                           "Neural Crest Cells","Migratory Neural Crest Cells"))


head(data@meta.data)
DefaultAssay(data) <- 'Spatial'
data <- NormalizeData(data,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(1,2,3)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)

