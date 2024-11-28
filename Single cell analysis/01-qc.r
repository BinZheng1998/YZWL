#!/usr/bin/env Rscript

# 使用 commandArgs 手动解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否包含 -h 或 --help 参数，并显示帮助信息
if ("-h" %in% args || "--help" %in% args) {
  cat("
Usage: Rscript script.R --input <path> --sample-name <name> --min-cells <int> --min-features <int> --percent-mito <float> [--min-feature-rna <int> --max-feature-rna <int> --min-feature-rna-percent <float> --max-feature-rna-percent <float> --min-mad-multiplier <float> --max-mad-multiplier <float>]

Options:
  --input                    Directory containing the 10X Genomics data (required)
  --sample-name              Name of the sample/project (required)
  --min-cells                Minimum number of cells expressing a gene (required)
  --min-features             Minimum number of genes detected in a cell (required)
  --percent-mito             Maximum percentage of mitochondrial genes allowed for quality control (required)
  --min-feature-rna          Minimum number of features (genes) for 'count' mode
  --max-feature-rna          Maximum number of features (genes) for 'count' mode
  --min-feature-rna-percent  Minimum feature percentage for 'percent' mode (0-1)
  --max-feature-rna-percent  Maximum feature percentage for 'percent' mode (0-1)
  --min-mad-multiplier       Minimum MAD multiplier for filtering
  --max-mad-multiplier       Maximum MAD multiplier for filtering

Example:
  Rscript script.R --input /path/to/data --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --min-mad-multiplier 2 --max-mad-multiplier 3
  ")
  quit(save = "no")
}

# 加载必要的库
library(Seurat)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(patchwork)

# 解析命令行参数
input_path <- args[which(args == "--input") + 1]
sample_name <- args[which(args == "--sample-name") + 1]
min_cells <- as.numeric(args[which(args == "--min-cells") + 1])
min_features <- as.numeric(args[which(args == "--min-features") + 1])
percent_mito <- as.numeric(args[which(args == "--percent-mito") + 1])

# 检查必要参数是否存在
if (is.null(input_path) || input_path == "") stop("Error: --input is missing or empty")
if (is.null(sample_name) || sample_name == "") stop("Error: --sample-name is missing or empty")

# 根据传入的参数设置过滤
min_feature_rna <- ifelse("--min-feature-rna" %in% args, as.numeric(args[which(args == "--min-feature-rna") + 1]), NA)
max_feature_rna <- ifelse("--max-feature-rna" %in% args, as.numeric(args[which(args == "--max-feature-rna") + 1]), NA)
min_feature_rna_percent <- ifelse("--min-feature-rna-percent" %in% args, as.numeric(args[which(args == "--min-feature-rna-percent") + 1]), NA)
max_feature_rna_percent <- ifelse("--max-feature-rna-percent" %in% args, as.numeric(args[which(args == "--max-feature-rna-percent") + 1]), NA)
min_mad_multiplier <- ifelse("--min-mad-multiplier" %in% args, as.numeric(args[which(args == "--min-mad-multiplier") + 1]), NA)
max_mad_multiplier <- ifelse("--max-mad-multiplier" %in% args, as.numeric(args[which(args == "--max-mad-multiplier") + 1]), NA)

# 打印参数确认
cat("Input path: ", input_path, "\n")
cat("Sample name: ", sample_name, "\n")

# 设置随机种子
set.seed(123)

# 判断输入的是数据目录还是 RDS 文件路径
if (grepl("\\.rds$", input_path)) {
  # 加载 RDS 文件
  seurat_obj <- readRDS(input_path)
  cat("Loaded Seurat object from RDS file.\n")
} else {
  # 加载 10X 数据
  counts <- Read10X(data.dir = input_path, gene.column = 1)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = min_cells, min.features = min_features)
  cat("Loaded 10X data from directory.\n")
}

# 添加线粒体基因百分比
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^gene-J6367")

# 创建输出的文件夹
output_dir <- paste0(sample_name, "_results")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  cat("Output directory created:", output_dir, "\n")
}
setwd(output_dir)

# 质控前的小提琴图
violin_p1 <- ggplot(seurat_obj@meta.data, aes(x =sample_name, y = nFeature_RNA)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle("nGene") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

violin_p2 <- ggplot(seurat_obj@meta.data, aes(x = sample_name, y = nCount_RNA)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle('nUMI') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

violin_p3 <- ggplot(seurat_obj@meta.data, aes(x = sample_name, y = percent.mt)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle('percent.mito') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

combined_violin <- violin_p1 | violin_p2 | violin_p3
ggsave(paste0(sample_name, "_violin_before_qc.png"), plot = combined_violin, width = 10, height = 5)

# 质控前的直方图
hist_p1 <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nGene", y = "Cell Counts") +
  theme_classic()

hist_p2 <- ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 5000, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nUMI", y = "Cell Counts") +
  theme_classic()

hist_p3 <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 5, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "percent.mito(%)", y = "Cell Counts") +
  scale_x_continuous(breaks = seq(0, max(seurat_obj@meta.data$percent.mt), by = 5)) +
  theme_classic()

combined_hist <- hist_p1 / hist_p2 / hist_p3
ggsave(paste0(sample_name, "_hist_before_qc.png"), plot = combined_hist, width = 10, height = 6)

# 质控前的散点图
calculate_correlation <- function(x, y) {
  test <- cor.test(x, y)
  cor_value <- round(test$estimate, 3)  # 相关系数
  p_value <- format(test$p.value, scientific = TRUE, digits = 3)
  return(list(cor_value = cor_value, p_value = p_value))
}

plot_with_correlation <- function(x, y, x_label, y_label, color = "blue") {
  cor_result <- calculate_correlation(x, y)
  cor_value <- cor_result$cor_value
  p_value <- cor_result$p_value

  ggplot(seurat_obj@meta.data, aes(x = x, y = y)) +
    geom_point(color = color, alpha = 0.9) +
    labs(x = x_label, y = y_label) +
    ggtitle(paste0("r=", cor_value, ", pvalue=", p_value)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}

scatter_plot1 <- plot_with_correlation(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$nFeature_RNA,
                                       "nUMI (nCount_RNA)", "nGene (nFeature_RNA)", "pink")
scatter_plot2 <- plot_with_correlation(seurat_obj@meta.data$nCount_RNA, seurat_obj@meta.data$percent.mt,
                                       "nUMI (nCount_RNA)", "percent.mt", "pink")
scatter_plot3 <- plot_with_correlation(seurat_obj@meta.data$nFeature_RNA, seurat_obj@meta.data$percent.mt,
                                       "nGene (nFeature_RNA)", "percent.mt", "pink")

combined_scatter <- scatter_plot1 | scatter_plot2 | scatter_plot3
ggsave(paste0(sample_name, "_scatter_before_qc.png"), plot = combined_scatter, width = 9, height = 3)

# 质控前的密度图
density_p1 <- ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nUMI", y = "Cell Counts") +
  theme_classic()

density_p2 <- ggplot(seurat_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nGene", y = "Cell Counts") +
  theme_classic()

density_p3 <- ggplot(seurat_obj@meta.data, aes(x = percent.mt)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "percent_mito", y = "Cell Counts") +
  theme_classic()

combined_density <- density_p1 | density_p2 | density_p3
ggsave(paste0(sample_name, "_density_before_qc.png"), plot = combined_density, width = 10, height = 5)

# 过滤数据
filter_conditions <- list()

if (!is.na(min_feature_rna)) {
  filter_conditions[["min_feature_rna"]] <- seurat_obj$nFeature_RNA > min_feature_rna
}
if (!is.na(max_feature_rna)) {
  filter_conditions[["max_feature_rna"]] <- seurat_obj$nFeature_RNA < max_feature_rna
}
if (!is.na(min_feature_rna_percent)) {
  filter_conditions[["min_feature_rna_percent"]] <- seurat_obj$nFeature_RNA / max(seurat_obj$nFeature_RNA) > min_feature_rna_percent
}
if (!is.na(max_feature_rna_percent)) {
  filter_conditions[["max_feature_rna_percent"]] <- seurat_obj$nFeature_RNA / max(seurat_obj$nFeature_RNA) < max_feature_rna_percent
}
if (!is.na(min_mad_multiplier)) {
  filter_conditions[["min_mad_multiplier"]] <- seurat_obj$nFeature_RNA > min_mad_multiplier
}
if (!is.na(max_mad_multiplier)) {
  filter_conditions[["max_mad_multiplier"]] <- seurat_obj$nFeature_RNA < max_mad_multiplier
}

# 合并过滤条件
combined_filter <- Reduce("&", filter_conditions)
seurat_obj_filtered <- subset(seurat_obj, subset = combined_filter & percent.mt < percent_mito)

# 保存过滤后的 Seurat 对象
saveRDS(seurat_obj_filtered, file = paste0(sample_name, "_filtered_seurat_obj.rds"))
cat("Filtered Seurat object saved: ", sample_name, "_filtered_seurat_obj.rds\n")

# 质控后的小提琴图
violin_p1 <- ggplot(seurat_obj_filtered@meta.data, aes(x = sample_name, y = nFeature_RNA)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle("nGene") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

violin_p2 <- ggplot(seurat_obj_filtered@meta.data, aes(x = sample_name, y = nCount_RNA)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle('nUMI') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

violin_p3 <- ggplot(seurat_obj_filtered@meta.data, aes(x = sample_name, y = percent.mt)) +
  geom_violin(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "", y = "") +
  ggtitle('percent.mito') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

combined_violin <- violin_p1 | violin_p2 | violin_p3
ggsave(paste0(sample_name, "_violin_after_qc.png"), plot = combined_violin, width = 10, height = 5)

# 过滤后的散点图
calculate_correlation <- function(x, y) {
  test <- cor.test(x, y)
  cor_value <- round(test$estimate, 3)  # 相关系数
   p_value <- format(test$p.value, scientific = TRUE, digits = 3)
  return(list(cor_value = cor_value, p_value = p_value))
}

plot_with_correlation <- function(x, y, x_label, y_label, color = "blue") {
  cor_result <- calculate_correlation(x, y)
  cor_value <- cor_result$cor_value
  p_value <- cor_result$p_value

  ggplot(seurat_obj_filtered@meta.data, aes(x = x, y = y)) +
    geom_point(color=color, alpha = 0.9) +
    labs(x = x_label, y = y_label) +
    ggtitle(paste0("r=", cor_value, ", pvalue=", p_value)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}

scatter_plot1 <- plot_with_correlation(seurat_obj_filtered@meta.data$nCount_RNA, seurat_obj_filtered@meta.data$nFeature_RNA,
                                       "nUMI (nCount_RNA)", "nGene (nFeature_RNA)", "pink")
scatter_plot2 <- plot_with_correlation(seurat_obj_filtered@meta.data$nCount_RNA, seurat_obj_filtered@meta.data$percent.mt,
                                       "nUMI (nCount_RNA)", "percent.mt", "pink")
scatter_plot3 <- plot_with_correlation(seurat_obj_filtered@meta.data$nFeature_RNA, seurat_obj_filtered@meta.data$percent.mt,
                                       "nGene (nFeature_RNA)", "percent.mt", "pink")
combined_scatter <- scatter_plot1 | scatter_plot2 | scatter_plot3
ggsave(paste0(sample_name, "_scatter_after_qc.png"), plot = combined_scatter, width = 9, height = 3)

# 过滤后的直方图
hist_p1 <- ggplot(seurat_obj_filtered@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nGene", y = "Cell Counts") +
  theme_classic()

hist_p2 <- ggplot(seurat_obj_filtered@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 5000, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nUMI", y = "Cell Counts") +
  theme_classic()

hist_p3 <- ggplot(seurat_obj_filtered@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 5, fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "percent.mito(%)", y = "Cell Counts") +
  scale_x_continuous(breaks = seq(0, max(seurat_obj@meta.data$percent.mt), by = 5)) +
  theme_classic()

combined_hist <- hist_p1 / hist_p2 / hist_p3
ggsave(paste0(sample_name, "_hist_after_qc.png"), plot = combined_hist, width = 10, height = 6)

# 过滤后的密度图
density_p1 <- ggplot(seurat_obj_filtered@meta.data, aes(x = nCount_RNA)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nUMI", y = "Cell Counts") +
  theme_classic()

density_p2 <- ggplot(seurat_obj_filtered@meta.data, aes(x = nFeature_RNA)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "nGene", y = "Cell Counts") +
  theme_classic()

density_p3 <- ggplot(seurat_obj_filtered@meta.data, aes(x = percent.mt)) +
  geom_density(fill = "pink", color = "black", alpha = 0.6) +
  labs(x = "percent_mito", y = "Cell Counts") +
  theme_classic()

combined_plot <- density_p1|density_p2|density_p3
ggsave(paste0(sample_name, "_density_after_qc.png"), plot = combined_plot, width = 10, height = 5)

cat("QC process and plotting completed.\n")
