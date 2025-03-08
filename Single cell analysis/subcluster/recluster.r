#!/usr/bin/env Rscript

# 加载所需包
suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(harmony)
    library(SeuratDisk)
    library(dplyr)
    library(patchwork)
    library(data.table)
    library(clustree)
    library(future)
    library(argparse)
})

# 定义命令行参数解析
parser <- ArgumentParser(description = "Process single-cell transcriptomics data with batch correction, clustering, and differential expression analysis.")
parser$add_argument("--input-rds", required = TRUE, 
                    help = "Path to input Seurat RDS file (required)")
parser$add_argument("--output-dir", default = "./", 
                    help = "Output directory [default: current directory './']")
parser$add_argument("--name", default = "single_cell_analysis", 
                    help = "Prefix for output files [default: 'single_cell_analysis']")
parser$add_argument("--cluster-id", default = NULL, 
                    help = "Comma-separated cluster IDs to subset (e.g., '1,2,3'), optional")
parser$add_argument("--batch", default = "replicate", 
                    help = "Batch variable for correction [default: 'replicate']")
parser$add_argument("--method", default = "integrate", choices = c("integrate", "harmony", "none"), 
                    help = "Batch correction method: 'integrate' (CCA), 'harmony', or 'none' [default: 'integrate']")
parser$add_argument("--dims", default = "1:15", 
                    help = "PCA dimensions for UMAP/clustering, R expression (e.g., '1:20') [default: '1:15']")
parser$add_argument("--resolution", default = 0.3, type = "double", 
                    help = "Clustering resolution, numeric value [default: 0.3]")

# 解析参数
args <- parser$parse_args()

# 将 dims 参数转换为整数向量
dims <- eval(parse(text = args$dims))

# 定义整合和聚类函数
run_Integrate <- function(obj, batch, dims, resolution) {
    obj_lst <- SplitObject(obj, split.by = batch)
    for (i in 1:length(obj_lst)) {
        DefaultAssay(obj_lst[[i]]) <- "RNA"
        obj_lst[[i]] <- NormalizeData(obj_lst[[i]], verbose = FALSE)
        obj_lst[[i]] <- FindVariableFeatures(obj_lst[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
    anchors <- FindIntegrationAnchors(object.list = obj_lst)
    combined <- IntegrateData(anchorset = anchors)
    combined <- ScaleData(combined, verbose = FALSE)
    combined <- RunPCA(combined, verbose = FALSE)
    combined <- RunUMAP(combined, dims = 1:dims, verbose = FALSE)
    combined <- FindNeighbors(combined, dims = 1:dims, verbose = FALSE)
    combined <- FindClusters(combined, resolution = resolution, verbose = FALSE)
    return(combined)
}

run_harmony <- function(obj, batch, dims, resolution) {
    obj <- SCTransform(obj, assay = "RNA", verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunHarmony(obj, group.by.vars = batch)
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:dims, verbose = FALSE)
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
    return(obj)
}

run_cluster <- function(obj, dims, resolution) {
    obj <- SCTransform(obj, assay = "RNA", verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:dims, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:dims, verbose = FALSE)
    obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
    return(obj)
}

# 主函数
process_single_cell <- function(input_rds, output_dir, name, cluster_id = NULL, batch, method, dims, resolution) {
    # 创建输出目录
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    # 读取数据
    obj <- readRDS(input_rds)
    
    # 提取子集（支持多个 cluster ID）
    if (!is.null(cluster_id)) {
        # 将 cluster_id 分割为向量
        cluster_ids <- unlist(strsplit(cluster_id, ","))
        # 转换为字符型并检查有效性
        obj@meta.data$seurat_clusters <- as.character(obj@meta.data$seurat_clusters)
        valid_clusters <- unique(obj@meta.data$seurat_clusters)
        invalid_ids <- setdiff(cluster_ids, valid_clusters)
        if (length(invalid_ids) > 0) {
            warning("The following cluster IDs are not found in the data: ", paste(invalid_ids, collapse = ", "))
        }
        # 提取匹配的子集
        obj <- subset(obj, seurat_clusters %in% cluster_ids)
    }

    # 根据方法进行批次校正和聚类
    if (method == "integrate") {
        obj <- run_Integrate(obj, batch, dims, resolution)
    } else if (method == "harmony") {
        obj <- run_harmony(obj, batch, dims, resolution)
    } else {
        obj <- run_cluster(obj, dims, resolution)
    }

    # UMAP 可视化
    p1 <- DimPlot(obj, group.by = "seurat_clusters", reduction = "umap", raster = FALSE, label = TRUE, repel = TRUE) + ggtitle("Clusters")
    p2 <- DimPlot(obj, group.by = "sample", reduction = "umap", raster = FALSE, label = TRUE, repel = TRUE) + ggtitle("Samples")
    png(file.path(output_dir, paste0(name, ".umap.png")), width = 12*300, height = 6*300, res = 300)
    print(p1 + p2)
    dev.off()

    # 聚类树
    sce_tmp <- FindClusters(obj, resolution = seq(0, 1.0, 0.1), verbose = FALSE)
    prefix <- ifelse(method == "integrate", "integrated_snn_res.", "SCT_snn_res.")
    p0 <- clustree(sce_tmp@meta.data, prefix = prefix) + labs(title = "Clustering Tree")
    png(file.path(output_dir, paste0(name, ".clustree.png")), width = 7*300, height = 12*300, res = 300)
    print(p0)
    dev.off()

    # 差异基因分析
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    plan(multisession, workers = 4)
    all.markers <- FindAllMarkers(obj, test.use = "wilcox", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    plan(sequential)
    write.csv(all.markers, file.path(output_dir, paste0(name, ".allmarkers.csv")))
    topn <- all.markers %>% filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
    write.csv(topn, file.path(output_dir, paste0(name, ".top50markers.csv")))

    # 保存处理后的对象
    saveRDS(obj, file.path(output_dir, paste0(name, ".processed.rds")))
}

# 执行函数
process_single_cell(
    input_rds = args$input_rds,
    output_dir = args$output_dir,
    name = args$name,
    cluster_id = args$cluster_id,
    batch = args$batch,
    method = args$method,
    dims = dims,
    resolution = args$resolution
)
