setwd('~/project/03_3D_spatial/02_result/202510-res/')
library(monocle3)
library(Seurat)
library(tidyverse)
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
options(future.globals.maxSize = 100* 1e10)
df <- readRDS("../../04_cell_annotation/new_data_20250827/E3.5_foregut.rds")

left_cell <- read.csv('E35_leftLung_cellID.csv')
left_cell$lung_type <- 'Left'
right_cell <- read.csv('E35_rightLung_cellID.csv')
right_cell$lung_type <- 'Right'

Lung_type <- rbind(left_cell,right_cell)
Lung_type <- Lung_type[,c(2,3)]

match_index <- match(rownames(df@meta.data), Lung_type$cell_barcode)
df@meta.data$lung_type <- Lung_type$lung_type[match_index]

type_vector <- setNames(Lung_type$lung_type, Lung_type$cell_barcode)
df <- AddMetaData(
  object = df,
  metadata = type_vector,
  col.name = "lung_type"
)
head(df@meta.data)

df <- UpdateSeuratObject(df)
df1 <- subset(df,subset = lung_type %in% c('Right','Left'))

all_markers <- FindAllMarkers(df1,test.use = 'DESeq2',group.by = 'lung_type',min.pct = 0.25,logfc.threshold = 0.25)
