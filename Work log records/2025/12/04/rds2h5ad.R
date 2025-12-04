setwd('~/project/03_3D_spatial/02_result/251203_SCENIC/')
library(Seurat)
library(SeuratDisk)
df1 <- readRDS('../../00_data/lung/E25_lung_bronchi.rds')
df1$celltype <- 'E25 Lung Epithliem'
df1$Stage <- 'E25'
df2 <- readRDS('../../00_data/lung/E35_lung_bronchi.rds')
df2$celltype <- 'E35 Lung Epithliem'
df2$Stage <- 'E35'
df3 <- readRDS('../../00_data/lung/E45_lung_bronchi.rds')
df3$celltype <- 'E45 Lung Epithliem'
df3$Stage <- 'E45'
df4 <- readRDS('../../00_data/lung/E25_lung_mesenchymal.rds')
df4$celltype <- 'E25 Lung Mesenchymal'
df4$Stage <- 'E25'
df5 <- readRDS('../../00_data/lung/E35_lung_mesenchymal_251129.rds')
df5$celltype <- 'E35 Lung Mesenchymal'
df5$Stage <- 'E35'
df6 <- readRDS('../../00_data/lung/E45_lung_mesenchymal_251129.rds')
df6$celltype <- 'E45 Lung Mesenchymal'
df6$Stage <- 'E45'
sample_ids <- c(
  "E25_Epi", "E35_Epi", "E45_Epi",
  "E25_Mes", "E35_Mes", "E45_Mes"
)
sc_combined <- merge(
  x = df1, 
  y = list(df2, df3, df4, df5, df6), 
  add.cell.ids = sample_ids, 
  project = "Chicken_Lung_Dev" 
)
obj <- sc_combined
obj <- JoinLayers(obj)
obj[["Spatial"]] <- as(object = obj[["Spatial"]], Class = "Assay")
DefaultAssay(obj) <- "Spatial"
obj <- NormalizeData(obj)
tmp <- DietSeurat(
  obj,
  assays = "Spatial",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE # SCENIC 不需要 misc
)

tmp@images <- list()
SaveH5Seurat(tmp, filename = "~/project/03_3D_spatial/00_data/lung/lung_251203.h5seurat", overwrite = TRUE)
Convert(
  "~/project/03_3D_spatial/00_data/lung/lung_251203.h5seurat", 
  assay = "Spatial", 
  dest = "h5ad", 
  overwrite = TRUE
)
