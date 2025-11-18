setwd('~/project/03_3D_spatial/02_result/251118_E35_E45_lung_mesenchymal_monocle3/')
library(Seurat)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(dplyr)
library(patchwork)
library(stringr)
library(rjson)
library(data.table)
library(harmony)
library(randomcoloR)

df <- readRDS("../251115_E35_lung_Mesenchymal_monocle3/E35_lung_mesenchymal.rds")
df1 <- readRDS("../251115_E45_lung_Mesenchymal_monocle3/E45_lung_mesenchymal_251116.rds")

####method1
table(df1$type)
df1$Stage[df1$type == "E45_single_slice"] <- 'E45_single_slice'
DimPlot(df1)
obj <- merge(df,df1)
table(obj$Stage)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj$type)  ## 对时期间和单张芯片去批次
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj,group.by = 'Stage',label = T,label.size = 8)

####method2
head(df)
df$type <- 'E35'
DimPlot(df,label = T)
df$type[df$Spatial_snn_res.1 == '5'] <- "E35_single_slice"

df$Stage <- 'E35'
df1$Stage <- 'E45'
obj <- merge(df,df1)
table(obj$Stage)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj$Stage)  ## 对时期间和单张芯片去批次
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- IntegrateLayers(object = obj, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj,label = T,label.size = 8)
DimPlot(obj,group.by = 'type',label = T,label.size = 8)
table(obj$type)

all_markers <- FindAllMarkers(object = obj,group.by = 'Spatial_snn_res.1',
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.5,
                              test.use = 'wilcox',
                              slot = 'data'
)
top50_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

####method3
table(obj$sample)
sample_counts <- table(obj$sample)
samples_to_keep <- names(sample_counts)[sample_counts > 5]
obj1 <- subset(obj, subset = sample %in% samples_to_keep)
obj1[["Spatial"]] <- JoinLayers(obj1[["Spatial"]])
obj1[["Spatial"]] <- split(obj1[["Spatial"]], f = obj1$sample)  
obj1 <- NormalizeData(obj1)
obj1 <- FindVariableFeatures(obj1)
obj1 <- ScaleData(obj1)
obj1 <- RunPCA(obj1)
obj1 <- FindNeighbors(obj1, dims = 1:30, reduction = "pca")
obj1 <- IntegrateLayers(object = obj1, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",
                       verbose = FALSE)
obj1[["Spatial"]] <- JoinLayers(obj1[["Spatial"]])
obj1 <- FindNeighbors(obj1, reduction = "integrated.harmony", dims = 1:30)
obj1 <- FindClusters(obj1, resolution = 1)
obj1 <- RunUMAP(obj1, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj1,label = T,label.size = 8)
DimPlot(obj1,group.by = 'sample',label = T,label.size = 8)
