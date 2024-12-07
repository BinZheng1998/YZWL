#!/usr/bin/env Rscript

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 帮助信息函数
show_help <- function() {
  cat("
Usage: Rscript test.r --input <input_file> --cont <contamination_threshold> --out <output_file> [--help|-h]

Arguments:
  --input                  Path to the input MTX data  (Seurat object)
  --cont                   Contamination threshold for filtering cells (default: 0.5)
  --out                    Path to the output RDS file (filtered Seurat object)

Options:
  --help, -h               Show this help message and exit

Example:
  Rscript test.r --input sample.rds --cont 0.5 --out filtered_sample.rds
")
}

# 检查是否需要显示帮助信息
if ("--help" %in% args || "-h" %in% args) {
  show_help()
  quit(save = "no", status = 0)
}

library(Seurat)
library(decontX)
library(Matrix)

# 提取参数值
get_arg <- function(flag, default = NULL) {
  pos <- which(args == flag)
  if (length(pos) > 0 && pos < length(args)) {
    return(args[pos + 1])
  }
  return(default)
}

# 获取命令行参数
input_file <- get_arg("--input")
cont_threshold <- as.numeric(get_arg("--cont", default = 0.5))  # 默认值 0.5
output_file <- get_arg("--out")

# 检查参数是否足够
if (is.null(input_file) || is.null(output_file)) {
  cat("Error: Missing required arguments.\n")
  show_help()
  quit(save = "no", status = 1)
}

# 打印参数信息
cat("Input file: ", input_file, "\n")
cat("Contamination threshold: ", cont_threshold, "\n")
cat("Output file: ", output_file, "\n")

# Step 1: 读取和合并数据
cat("Loading data...\n")
counts <- Read10X(data.dir = input_file, gene.column = 1)
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
seurat_object2 <- JoinLayers(seurat_obj)

# Step 2: 提取计数矩阵
cat("Extracting counts matrix...\n")
counts_matrix <- GetAssayData(object = seurat_object2,assay="RNA", layer = "counts")
if (!inherits(counts_matrix, "dgCMatrix")) {
  counts_matrix <- as(counts_matrix, "dgCMatrix")
}

# Step 3: 去污染处理
cat("Running decontX...\n")
decontx_results <- decontX(counts_matrix)

# Step 4: 添加污染信息
cat("Adding contamination results...\n")
seurat_object2$contamination <- decontx_results$contamination

# Step 5: 数据过滤
cat("Filtering cells based on contamination threshold...\n")
filtered_seurat_object <- subset(seurat_object2, subset = contamination < cont_threshold)

# Step 6: 保存结果
cat("Saving filtered data to output file...\n")
saveRDS(filtered_seurat_object, output_file)

cat("Process completed successfully!\n")
