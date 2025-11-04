setwd('~/project/03_3D_spatial/02_result/251104_lung_progenitor_cells_cellcycle/')
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(patchwork)
library(scplotter)

set.seed(12345)
options(future.globals.maxSize = 100* 1e10)

obj_sub <- readRDS("../251101_lung_hdWGCNA/E35_45_lung_progenitor_cells.rds")

s.genes <- read.table('~/project/04_gonad_single_cell/new_result/script/chicken_s.gene.txt')
s.genes$V1 <- gsub('gene-','',s.genes$V1)
g2m.genes <- read.table('~/project/04_gonad_single_cell/new_result/script/chicken_g2m.gene.txt')
g2m.genes$V1 <- gsub('gene-','',g2m.genes$V1)

res <- CellCycleScoring(obj_sub,s.features = s.genes$V1,g2m.features = g2m.genes$V1,set.ident = T)
head(res)

pdf('E35_45_lung_progenitor_cells_cellcycle.pdf')
CellDimPlot(res,group_by = 'Phase',theme= ggplot2::theme_classic, theme_args = list(base_size =16),palette ="seurat")
dev.off()

pdf('E35_45_lung_progenitor_cells_cellcycle_G2Mscore.pdf')
FeatureStatPlot(res, plot_type ="dim", 
                features ="G2M.Score", 
                reduction ="umap",
                theme= ggplot2::theme_classic,
                theme_args = list(base_size =16))
dev.off()

pdf('E35_45_lung_progenitor_cells_cellcycle_Sscore.pdf')
FeatureStatPlot(res, plot_type ="dim", 
                features ="S.Score", 
                reduction ="umap",
                theme= ggplot2::theme_classic,
                theme_args = list(base_size =16))
dev.off()
