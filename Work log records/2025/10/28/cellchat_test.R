setwd('~/project/03_3D_spatial/02_result/2510_cellchat/')
library(Seurat)
library(CellChat)
library(SeuratData)#加载seurat数据集
library(tidyverse)
packageVersion("CellChat")
load('CellChatDB.rda')

#InstallData("pbmc3k")
#data("pbmc3k")
data@meta.data$cell <- data@meta.data$Spatial_snn_res.0.3
sce <- UpdateSeuratObject(data)
table(sce$cell)
colnames(sce@meta.data)
dim(sce)
# 去掉没有注释信息的细胞
sce <- sce[ ,which(!is.na(sce@meta.data$cell))]
sce <- sce %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData

sce
table(Idents(sce))
Idents(sce) <-"cell"


## 1.输入数据
# For the gene expression data matrix, genes should be in rows with rownames and cells in columns with colnames.
# Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required
# as input for CellChat analysis
data.input = sce@assays$Spatial@data# normalized data matrix
data.input[1:4,1:4]
meta = sce@meta.data[,"cell",drop=F]# a dataframe with rownames containing cell mata data
colnames(meta) <-"labels"
meta$samples <- 'sample1'
meta$samples <- factor(meta$samples)
meta$labels <- paste0('c',meta$labels)
head(meta)
table(meta)

## 2.创建对象
cellchat <- createCellChat(object = data.input, meta = meta, group.by ="labels")
cellchat
levels(cellchat@idents)# show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))# number of cells in each cell group
groupSize


## 3.数据库
#CellChatDB <- CellChatDB.human# use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat


## 4.鉴定亚群高表达基因
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat)# This step is necessary even if using the whole database
future::plan("multisession", workers = 20)# do parallel
# CellChat identifies over-expressed ligands or receptors in one cell group
cellchat <- updateCellChat(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat

## 5.计算probability
cellchat <- computeCommunProb(cellchat,type="triMean")
#> triMean is used for calculating the average gene expression per cell group.

## 6.通路水平的通讯
cellchat <- computeCommunProbPathway(cellchat)

## 7.计算汇总的通讯网络
cellchat <- aggregateNet(cellchat)

## 8.提取细胞通讯结果
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors
# 默认 thresh ：threshold of the p-value for determining significant interaction
df.net <- subsetCommunication(cellchat, thresh = 0.05)
head(df.net)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show <-"WNT"
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap ="Reds")

