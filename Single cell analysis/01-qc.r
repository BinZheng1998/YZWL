#!/usr/bin/env Rscript

# 使用 commandArgs 手动解析命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否包含 -h 或 --help 参数，并显示帮助信息
if ("-h" %in% args || "--help" %in% args) {
  cat("
Usage: Rscript script.R --sample-dir <directory> --sample-name <name> --min-cells <int> --min-features <int> --percent-mito <float> --filter-mode <mode> [--min-feature-rna <int> --max-feature-rna <int> --min-feature-rna-percent <float> --max-feature-rna-percent <float> --mad-multiplier <float>]

Options:
  --sample-dir               Directory containing the 10X Genomics data (required)
  --sample-name              Name of the sample/project (required)
  --min-cells                Minimum number of cells expressing a gene (required)
  --min-features             Minimum number of genes detected in a cell (required)
  --percent-mito             Maximum percentage of mitochondrial genes allowed for quality control (required)
  --filter-mode              Filtering mode: 'count', 'percent', or 'mad' (required)
  --min-feature-rna          Minimum number of features (genes) for 'count' mode
  --max-feature-rna          Maximum number of features (genes) for 'count' mode
  --min-feature-rna-percent  Minimum feature percentage for 'percent' mode (0-1)
  --max-feature-rna-percent  Maximum feature percentage for 'percent' mode (0-1)
  --mad-multiplier           Multiplier for MAD-based filtering thresholds for 'mad' mode

Example:
  Rscript script.R --sample-dir /path/to/data --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --filter-mode mad --mad-multiplier 2  Rscript script.R --sample-dir /path/to/data --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --filter-mode count --min-feature-rna 200 --max-feature-rna 5000
  Rscript script.R --sample-dir /path/to/data --sample-name sample1 --min-cells 3 --min-features 200 --percent-mito 5 --filter-mode percent --min-feature-rna-percent 0.1 --max-feature-rna-percent 0.9
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
data_dir <- args[which(args == "--sample-dir") + 1]
sample_name <- args[which(args == "--sample-name") + 1]
min_cells <- as.numeric(args[which(args == "--min-cells") + 1])
min_features <- as.numeric(args[which(args == "--min-features") + 1])
percent_mito <- as.numeric(args[which(args == "--percent-mito") + 1])
filter_mode <- args[which(args == "--filter-mode") + 1]

# 检查必要参数是否存在
if (is.null(data_dir) || data_dir == "") stop("Error: --sample-dir is missing or empty")
if (is.null(sample_name) || sample_name == "") stop("Error: --sample-name is missing or empty")
if (is.null(filter_mode)) stop("Error: --filter-mode is missing or empty")
if (!filter_mode %in% c("count", "percent", "mad")) stop("Error: --filter-mode must be 'count', 'percent', or 'mad'")

# 根据模式解析参数
if (filter_mode == "count") {
  min_feature_rna <- as.numeric(args[which(args == "--min-feature-rna") + 1])
  max_feature_rna <- as.numeric(args[which(args == "--max-feature-rna") + 1])
} else if (filter_mode == "percent") {
  min_feature_rna_percent <- as.numeric(args[which(args == "--min-feature-rna-percent") + 1])
  max_feature_rna_percent <- as.numeric(args[which(args == "--max-feature-rna-percent") + 1])
} else if (filter_mode == "mad") {
  mad_multiplier <- as.numeric(args[which(args == "--mad-multiplier") + 1])
}

# 打印参数确认
cat("Sample directory: ", data_dir, "\n")
cat("Sample name: ", sample_name, "\n")
cat("Filter mode: ", filter_mode, "\n")

# 设置随机种子
set.seed(123)

# 使用 Read10X 直接读取数据
counts <- Read10X(data.dir = data_dir, gene.column = 1)

# 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = min_cells, min.features = min_features)

# 添加线粒体基因百分比
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^gene-J6367")
#seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^gene-RP[SL]")

#创建输出的文件夹
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


# 根据模式进行过滤
if (filter_mode == "count") {
  cat("Applying count-based filtering...\n")
  seurat_obj <- subset(seurat_obj, subset = 
                         nFeature_RNA > min_feature_rna & nFeature_RNA < max_feature_rna & 
                         percent.mt < percent_mito)
} else if (filter_mode == "percent") {
  cat("Applying percent-based filtering...\n")
  feature_quantiles <- quantile(seurat_obj$nFeature_RNA, probs = c(min_feature_rna_percent, max_feature_rna_percent))
  seurat_obj <- subset(seurat_obj, subset = 
                         nFeature_RNA > feature_quantiles[1] & nFeature_RNA < feature_quantiles[2] & 
                         percent.mt < percent_mito)
} else if (filter_mode == "mad") {
  cat("Applying MAD-based filtering...\n")
  nFeature_median <- median(seurat_obj$nFeature_RNA)
  nFeature_mad <- mad(seurat_obj$nFeature_RNA)
  seurat_obj <- subset(seurat_obj, subset = 
                         nFeature_RNA > (nFeature_median - mad_multiplier * nFeature_mad) & 
                         nFeature_RNA < (nFeature_median + mad_multiplier * nFeature_mad) & 
                         percent.mt < percent_mito)
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


saveRDS(seurat_obj, file = paste0(sample_name, "_after_qc.rds"))

cat("Data processing complete. Filtered data saved to: ", paste0(sample_name, "_after_qc.rds"), "\n")
