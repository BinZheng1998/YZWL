
setwd('~/project/03_3D_spatial/02_result/260121_cellchat/')
library(Seurat)
library(CellChat)
library(SeuratData)
library(tidyverse)
packageVersion("CellChat")
options(future.globals.maxSize = 200 * 1024^3)
load('~/project/03_3D_spatial/00_data/database/cellchat/CellChatDB.chicken.rda')

df <- readRDS('../../00_data/lung/E45_lung_epithelium_251220.rds')
head(df)
df1 <- readRDS('../../00_data/lung/E45_lung_mesenchymal_251220.rds')
df1$lung_celltype[df1$lung_celltype == "Distal Mesenchymal"] <- "Distal Mesenchyme"
df1$lung_celltype[df1$lung_celltype == "Proximal Mesenchymal"] <- "Proximal Mesenchyme"
df1$celltype <- df1$lung_celltype
head(df1)

sce <- merge(x = df,y = df1, 
                       add.cell.ids = c("Epithelium", "Mesenchyme"), 
                       project = "E45_Lung")

sce <- UpdateSeuratObject(sce)
sce@meta.data$cell <- sce@meta.data$celltype

sce <- sce %>%
  NormalizeData %>%
  FindVariableFeatures %>%
  ScaleData

Idents(sce) <-"cell"
sce
sce <- JoinLayers(sce,assay = 'Spatial')

CellChatDB.use <- subsetDB(CellChatDB.chicken)
#cellchat@DB <- CellChatDB.use


scale.ratio <- 4.125
spot.radius <- scale.ratio / 1
spatial.factors <- data.frame(
  ratio = scale.ratio, 
  tol = 2.1
)

createCellChat_3D <-function(object, meta = NULL, group.by = NULL,
                              datatype = c("RNA", "spatial"), 
                              coordinates = NULL, spatial.factors = NULL,
                              assay = NULL, do.sparse = T) {
  
  datatype <- match.arg(datatype)
  
  # --- 1. 处理输入数据类型 (Seurat / Matrix / SCE) ---
  if (inherits(x = object, what = c("matrix", "Matrix", "dgCMatrix", "dgRMatrix","CsparseMatrix"))) {
    message("Create a CellChat object from a data matrix")
    data <- object
    if (is.null(group.by)) { group.by <- "labels" }
  }
  
  if (is(object,"Seurat")) {
    message("Create a CellChat object from a Seurat object")
    if (is.null(assay)) {
      assay = Seurat::DefaultAssay(object)
    }
    # 兼容 Seurat v3/v4 和 v5
    if (packageVersion("Seurat") < "5.0.0") {
      data <- object[[assay]]@data
    } else {
      data <- object[[assay]]$data
    }
    
    if (min(data) < 0) { stop("The data matrix contains negative values.") }
    
    if (is.null(meta)) {
      meta <- object@meta.data
      meta$ident <- Seurat::Idents(object)
    }
    if (is.null(group.by)) { group.by <- "ident" }
    
    # 自动获取坐标 (如果是 Seurat 且没提供坐标)
    if (datatype %in% c("spatial") && is.null(coordinates)) {
      coordinates <- Seurat::GetTissueCoordinates(object, scale = NULL, cols = c("imagerow", "imagecol"))
    }
  }
  
  # --- 2. 处理 Metadata ---
  if (!is.null(meta)) {
    if (!is.data.frame(meta)) { meta <- as.data.frame(x = meta) }
    # 强制行名对齐
    if (!identical(rownames(meta), colnames(data))) {
      warning("Barcodes in meta do not match data columns. Forcing assignment.")
      rownames(meta) <- colnames(data)
    }
  } else {
    meta <- data.frame()
  }
  
  # --- 3. 核心修改：支持 3D 坐标 ---
  if (datatype %in% c("spatial")) {
    
    # 【修改点开始】：允许 2 列或 3 列输入
    cols <- ncol(coordinates)
    if (cols == 2) {
      colnames(coordinates) <- c("x_cent","y_cent")
    } else if (cols == 3) {
      colnames(coordinates) <- c("x_cent","y_cent", "z_cent")
      message(">>> 检测到 3D 坐标 (x, y, z)。将启用 3D 距离计算。注意：部分绘图函数可能不兼容。<<<")
    } else {
      stop("Please check the input 'coordinates'. It must be a 2-column (XY) or 3-column (XYZ) matrix.")
    }
    # 【修改点结束】
    
    if (is.null(spatial.factors) | !("ratio" %in% names(spatial.factors)) | !("tol" %in% names(spatial.factors))) {
      stop("spatial.factors with colnames `ratio` and `tol` should be provided!")
    } else {
      coordinates <- as.matrix(coordinates)
      images = list("coordinates" = coordinates, "spatial.factors" = spatial.factors)
    }
    message("Create a CellChat object from spatial transcriptomics data...")
  } else {
    images <- list()
  }
  
  # --- 4. 创建对象 ---
  # 注意：这里需要直接调用 new()，因为我们修改了逻辑
  object_cc <- new(Class = "CellChat",
                   data = data,
                   images = images,
                   meta = meta)
  
  # --- 5. 设置 Group 和 Sample ---
  if (!is.null(meta) & nrow(meta) > 0) {
    if (!("samples" %in% colnames(meta))) {
      meta$samples <- "sample1"
    }
    meta$samples <- factor(meta$samples)
    object_cc@meta <- meta
    
    if (!(group.by %in% colnames(meta))) {
      stop("The 'group.by' column is missing in metadata.")
    }
    object_cc <- setIdent(object_cc, ident.use = group.by)
  }
  
  object_cc@options$mode <- "single"
  object_cc@options$datatype <- datatype
  
  return(object_cc)
}



process_one_sample <- function(seu_obj, name_id, db_use, spatial_factors) {
  
  message(paste0(">>> Processing Sample: ", name_id, " <<<"))
  
  # 1. 读取数据
  data.input <- Seurat::GetAssayData(seu_obj, layer = "data", assay = "Spatial")
  meta <- seu_obj@meta.data
  
  # 2. 标签处理
  if("seurat_clusters3" %in% colnames(meta)){
    meta$labels <- meta$seurat_clusters3
  } else {
    meta$labels <- Idents(seu_obj) 
  }
  meta$samples <- name_id
  
  # 3. 智能坐标提取
  possible_combinations <- list(
    suffix_x = c("trans_x.x", "trans_y.x", "trans_z.x"), 
    suffix_y = c("trans_x.y", "trans_y.y", "trans_z.y"), 
    no_suffix = c("trans_x", "trans_y", "trans_z"),      
    standard  = c("x", "y", "z")                         
  )
  
  coord_cols <- NULL
  found_type <- ""
  
  for (type in names(possible_combinations)) {
    cols <- possible_combinations[[type]]
    if (all(cols %in% colnames(meta))) {
      coord_cols <- cols
      found_type <- type
      break 
    }
  }
  
  if (is.null(coord_cols)) {
    stop(paste0("❌ Error: 样本 ", name_id, " 找不到坐标列！"))
  } else {
    message(paste0("✅ Found coordinates type: ", found_type, " -> ", paste(coord_cols, collapse=", ")))
  }
  
  locs <- meta[, coord_cols]
  colnames(locs) <- c("x", "y", "z")
  
  # 4. 创建 CellChat (Stereo-seq Bin50)
  cellchat <- createCellChat_3D(
    object = data.input, 
    meta = meta, 
    group.by = "labels", 
    datatype = "spatial", 
    coordinates = as.matrix(locs), 
    spatial.factors = spatial.factors
  )
  cellchat@DB <- db_use
  
  # 5. 分析流程
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # -------------------------------------------------------------
  # 6. 核心参数修正
  # -------------------------------------------------------------
  
  # 检查 spatial.factors 是否包含 ratio，如果有则使用，否则为 NULL
  my_scale <- if(!is.null(spatial.factors) && "ratio" %in% names(spatial.factors)) spatial.factors$ratio else NULL
  
  message(paste0(">>> 计算通讯概率参数: Range = 250 um, Scale Factor = ", my_scale, " <<<"))
  
  cellchat <- computeCommunProb(
    cellchat, 
    type = "truncatedMean", 
    trim = 0.1, 
    
    distance.use = TRUE,    
    interaction.range = 250, 
    scale.distance = my_scale, 
    
    contact.dependent = TRUE, 
    contact.range = 100
  )
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}

cellchat.E4.5 <- process_one_sample(
  seu_obj = sce,       # 输入你的 Seurat 对象变量名
  name_id = "E4.5",        # 样本名称
  db_use = CellChatDB.use,
  spatial_factors = spatial.factors
)

saveRDS(cellchat.E4.5, file = "cellchat_E4.5_spatial_analysis2.rds")
