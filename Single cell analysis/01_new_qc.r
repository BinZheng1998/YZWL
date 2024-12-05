#!/usr/bin/env Rscript

# 使用 commandArgs 手动解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否包含 -h 或 --help 参数，并显示帮助信息
if ("-h" %in% args || "--help" %in% args) {
  cat("
Usage: Rscript script.R --input <path> --sample-name <name> --min-cells <int> --min-features <int> --percent-mito <float> --filter-mode <mode> [--min-feature-rna <int> --max-feature-rna <int> --min-feature-rna-percent <float> --max-feature-rna-percent <float> --min-mad-multiplier <float> --max-mad-multiplier <float>]

Options:
  --input                   Path to the 10X Genomics data directory or an RDS file (required)
  --sample-name             Name of the sample/project (required)
  --min-cells               Minimum number of cells expressing a gene (required)
  --min-features            Minimum number of genes detected in a cell (required)
  --percent-mito            Maximum percentage of mitochondrial genes allowed for quality control (required)
  --filter-mode             Filtering mode:
                              'mode1': Use --min-feature-rna and --max-feature-rna
                              'mode2': Use --min-feature-rna and --max-feature-rna-percent
                              'mode3': Use --min-feature-rna and --max-mad-multiplier
                              'mode4': Use --max-feature-rna and --min-mad-multiplier
                              'mode5': Use --min-feature-rna-percent and --max-feature-rna
                              'mode6': Use --min-mad-multiplier and --max-mad-multiplier
  --min-feature-rna         Minimum number of features (genes) for some modes
  --max-feature-rna         Maximum number of features (genes) for some modes
  --min-feature-rna-percent Minimum feature percentage for some modes (0-1)
  --max-feature-rna-percent Maximum feature percentage for some modes (0-1)
  --min-mad-multiplier      Minimum multiplier for MAD-based filtering for some modes
  --max-mad-multiplier      Maximum multiplier for MAD-based filtering for some modes

Examples:
  Rscript script.R --input /path/to/data --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --filter-mode mode1 --min-feature-rna 200 --max-feature-rna 5000
  Rscript script.R --input /path/to/data/sample1_after_qc.rds --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --filter-mode mode6 --min-mad-multiplier 1 --max-mad-multiplier 2
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
filter_mode <- args[which(args == "--filter-mode") + 1]
min_feature_rna <- as.numeric(args[which(args == "--min-feature-rna") + 1])
max_feature_rna <- as.numeric(args[which(args == "--max-feature-rna") + 1])
min_feature_rna_percent <- as.numeric(args[which(args == "--min-feature-rna-percent") + 1])
max_feature_rna_percent <- as.numeric(args[which(args == "--max-feature-rna-percent") + 1])
min_mad_multiplier <- as.numeric(args[which(args == "--min-mad-multiplier") + 1])
max_mad_multiplier <- as.numeric(args[which(args == "--max-mad-multiplier") + 1])

cat("Sample Path: ", input_path, "\n")
cat("Sample name: ", sample_name, "\n")
cat("Filter mode: ", filter_mode, "\n")

# 设置随机种子
set.seed(123)

# 根据输入类型读取数据
if (grepl("\\.rds$", input_path, ignore.case = TRUE)) {
  # 读取RDS文件
  seurat_obj <- readRDS(input_path)
} else {
  # 读取10X Genomics数据
  counts <- Read10X(data.dir = input_path)
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = min_cells, min.features = min_features)
}

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^gene-J6367")

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

#质控

if (filter_mode == "mode1") {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_feature_rna & nFeature_RNA <= max_feature_rna & percent.mt <= percent_mito)
} else if (filter_mode == "mode2") {
  max_rna <- quantile(seurat_obj$nFeature_RNA, probs = max_feature_rna_percent)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_feature_rna & nFeature_RNA <= max_rna & percent.mt <= percent_mito)
} else if (filter_mode == "mode3") {
  mad_val <- mad(seurat_obj$nFeature_RNA, constant = 1) * max_mad_multiplier
  median_val <- median(seurat_obj$nFeature_RNA)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_feature_rna & nFeature_RNA <= (median_val + mad_val) & percent.mt <= percent_mito)
} else if (filter_mode == "mode4") {
  mad_val <- mad(seurat_obj$nFeature_RNA, constant = 1) * min_mad_multiplier
  median_val := median(seurat_obj$nFeature_RNA)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= (median_val - mad_val) & nFeature_RNA <= max_feature_rna & percent.mt <= percent_mito)
} else if (filter_mode == "mode5") {
  min_rna <- quantile(seurat_obj$nFeature_RNA, probs = min_feature_rna_percent)
  max_rna := quantile(seurat_obj$nFeature_RNA, probs = max_feature_rna_percent)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= min_rna & nFeature_RNA <= max_rna & percent.mt <= percent_mito)
} else if (filter_mode == "mode6") {
  mad_min := mad(seurat_obj$nFeature_RNA, constant = 1) * min_mad_multiplier
  mad_max := mad(seurat_obj$nFeature_RNA, constant = 1) * max_mad_multiplier
  median_val := median(seurat_obj$nFeature_RNA)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA >= (median_val - mad_min) & nFeature_RNA <= (median_val + mad_max) & percent.mt <= percent_mito)
}

# 质控后的小提琴图
violin_p1 <- ggplot(seurat_obj@meta.data, aes(x = sample_name, y = nFeature_RNA)) +
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

  ggplot(seurat_obj@meta.data, aes(x = x, y = y)) +
    geom_point(color=color, alpha = 0.9) +
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
ggsave(paste0(sample_name, "_scatter_after_qc.png"), plot = combined_scatter, width = 9, height = 3)

# 过滤后的直方图
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
ggsave(paste0(sample_name, "_hist_after_qc.png"), plot = combined_hist, width = 10, height = 6)

# 过滤后的密度图
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

combined_plot <- density_p1|density_p2|density_p3
ggsave(paste0(sample_name, "_density_after_qc.png"), plot = combined_plot, width = 10, height = 5)

# 质控后保存结果
saveRDS(seurat_obj, file = paste0(sample_name, "_qc.rds"))
cat("Data processing complete. Filtered data saved to: ", paste0(sample_name, "_qc.rds"), "\n")
