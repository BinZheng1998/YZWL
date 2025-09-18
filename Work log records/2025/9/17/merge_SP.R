setwd(dir = '~/project/03_3D_spatial/04_cell_annotation/new_data_2509/merge/')
library(Seurat)
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
options(future.globals.maxSize = 1e10)
packageVersion("Seurat")
obj.list = list()
stage_list <- c('E2.5','E3.5','E4.5')
for (i in 1:length(stage_list)){
  obj.list[[i]] <- readRDS(paste0('../../new_data_20250827/',stage_list[i],'_foregut.rds'))
  obj.list[[i]]$Stage <- stage_list[i]
}

#不合并
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, verbose = F)
obj <- RunPCA(obj, verbose = F)
obj <- FindNeighbors(obj, reduction = 'pca', dims = 1:30)
obj <- FindClusters(obj, resolution = 1.0)
obj <- RunUMAP(obj, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
p1 <-DimPlot(obj, group.by='Stage', reduction = 'umap', raster=FALSE,label=F,repel =T)
p1
head(obj@meta.data)
p2 <-DimPlot(obj, group.by='celltype_gut', reduction = 'umap' ,raster=FALSE,label=F,repel =T,label.size=2)+
  theme(legend.position = 'none')
p2

#CCA合并
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj$Stage)  ## 对时期间去批次
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")

obj <- IntegrateLayers(object = obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:30)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, reduction = 'integrated.cca', dims = 1:30, return.model = T, verbose = F)

options(repr.plot.width =22,repr.plot.height =7)
DimPlot(obj, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')


#RPCA合并
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj$Stage)  ## 对时期间去批次
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")

obj <- IntegrateLayers(object = obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                       verbose = FALSE)
obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:30)
obj <- FindClusters(obj, resolution = 1)
obj <- RunUMAP(obj, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)

options(repr.plot.width =22,repr.plot.height =7)
DimPlot(obj, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')

#harmony合并
obj <- merge(obj.list[[1]],obj.list[2:length(obj.list)])
obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj$Stage)  ## 对时期间去批次
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
DimPlot(obj, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')






####Lung
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                'Lung Progenitor Cells'))
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj_sub <- ScaleData(obj_sub, verbose = F)
obj_sub <- RunPCA(obj_sub, verbose = F)
obj_sub <- FindNeighbors(obj_sub, reduction = 'pca', dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1.0)
obj_sub <- RunUMAP(obj_sub, reduction = 'pca', dims = 1:30, return.model = T, verbose = F)
p1 <-DimPlot(obj_sub, group.by='Stage', reduction = 'umap', raster=FALSE,label=F,repel =T)
p1
p2 <-DimPlot(obj_sub, group.by='celltype_gut', reduction = 'umap' ,raster=FALSE,label=F,repel =T,label.size=2)
  theme(legend.position = 'none')
p2


#Lung CCA
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                   'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = CCAIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.cca",
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.cca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.cca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))

#Lung RPCA
options(future.globals.maxSize = 100 * 1024^3) 
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                   'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                           new.reduction = "integrated.rpca",
                           k.weight = 50,k.anchor = 1,
                           verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))

#Lung harmony
obj_sub <- subset(obj,subset = celltype_gut %in% c('Lung Mesenchymal Progenitor Cells','Lung Bud1','Lung Bud2',
                                                   'Lung Secondary Bronchi','Lung Primary Bronchi',
                                                   'Lung Progenitor Cells'))
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)  ## 对时期间去批次
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = HarmonyIntegration, orig.reduction = "pca", 
                           new.reduction = "integrated.harmony",theta=0.8,lambda=2,sigma=0.2,
                           verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.harmony", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub, reduction = "umap", group.by = c("Stage", "celltype_gut"))+theme(legend.position = 'none')
DimPlot(obj_sub, reduction = "umap", group.by = c("celltype_gut"))
#marker_genes <- c('NKX2-1','MME','WNT4','NKX2-8','RIPPLY3','SYT10','MLLT3','ADCY8','TBX4')
#raw=FetchData(obj_sub,c(marker_genes,"celltype_gut"))
#harmony_res=FetchData(obj_sub,c(marker_genes,"celltype_gut"),layer = 'counts')
