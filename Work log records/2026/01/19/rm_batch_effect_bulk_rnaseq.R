# ==============================================================================
# Part 1: 全局准备 (Global Preparation) - 只运行一次
# ==============================================================================
setwd('~/project/04_rnaseq/pig/2.result/260119_rmbatcheffect/')
library(DESeq2)
library(sva)
library(limma)
library(ggplot2)
library(patchwork)

# --- 1.1 读取和清洗原始数据 ---
cat("正在读取大矩阵...\n")
df <- read.table('../../1.data/pig_total_merged_counts_matrix.txt', sep = '\t', header = T)
rownames(df) <- df$Geneid
df1 <- df[,-1]

# 过滤低质量样本 (Total Counts < 1M)
df2 <- df1[, colSums(df1, na.rm = T) >= 1000000]

# 过滤离群样本 (基于检出基因数)
detected_genes_gt1 <- colSums(df2 > 1)
Q1 <- quantile(detected_genes_gt1, 0.25)
Q3 <- quantile(detected_genes_gt1, 0.75)
IQR_val <- IQR(detected_genes_gt1)
whisker_lower <- max(0, Q1 - 1.5 * IQR_val) 
whisker_upper <- Q3 + 1.5 * IQR_val

# 仅保留合格样本
samples_to_keep <- names(detected_genes_gt1[detected_genes_gt1 >= whisker_lower])
df3 <- df2[, samples_to_keep]

# --- 1.2 整理 Metadata ---
metadata <- read.table('../../1.data/metadata/pig_metadata.txt', sep = '\t', header = T)
metadata1 <- metadata[, c(2, 3, 9, 19)]
colnames(metadata1) <- c('Sample', 'Project', 'BW', 'Tissue')

# 筛选需要的组织和分组
target_tissues <- c("Liver", "Muscle", "Hypothalamus", "Adipose", "Pituitary")
metadata2 <- metadata1[metadata1$Tissue %in% target_tissues, ]
metadata2 <- metadata2[metadata2$BW %in% c("High", "Low"), ]
metadata2$Group <- paste0(metadata2$BW, sep = '_', metadata2$Tissue)

# 定义 PCA 绘图函数 (保持不变，用于出图)
plot_pca <- function(mat, info, title) {
  pca <- prcomp(t(mat))
  df_pca <- as.data.frame(pca$x)
  info <- info[colnames(mat), ] 
  df_pca$Group <- info$Group
  df_pca$Project <- info$Project
  
  perc <- round(100 * pca$sdev^2 / sum(pca$sdev^2))
  
  ggplot(df_pca, aes(PC1, PC2, shape = Group, color = Project)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = title, 
         x = paste0("PC1 (", perc[1], "%)"), 
         y = paste0("PC2 (", perc[2], "%)")) +
    theme_bw()
}

# ==============================================================================
# Part 2: 循环处理每个器官 (Loop Execution)
# ==============================================================================

for (tissue in target_tissues) {
  
  message(paste0("\n========================================================"))
  message(paste0("正在处理器官: ", tissue))
  message(paste0("========================================================"))
  
  # --- 2.1 针对当前器官取子集 ---
  metadata_sub <- metadata2[metadata2$Tissue == tissue, ]
  valid_cols <- intersect(metadata_sub$Sample, colnames(df3))
  
  if(length(valid_cols) < 4) {
    message(paste("跳过:", tissue, "- 样本量太少 (<4)"))
    next
  }
  
  df_sub <- df3[, valid_cols]
  coldata_sub <- metadata_sub[metadata_sub$Sample %in% valid_cols, ]
  rownames(coldata_sub) <- coldata_sub$Sample
  df_sub <- df_sub[, rownames(coldata_sub)]
  
  # --- 2.2 基因过滤 (已修正 bug) ---
  # 动态调整过滤阈值
  #min_samples <- min(10, floor(ncol(df_sub)/2))
  keep <- rowSums(df_sub > 5) >= 10
  df_final <- df_sub[keep, ]
  
  message(paste("保留基因数:", nrow(df_final)))
  
  # 动态设置分组
  group_high <- paste0("High_", tissue)
  group_low  <- paste0("Low_", tissue)
  coldata_sub$Group <- factor(coldata_sub$Group, levels = c(group_low, group_high))
  contrast_vec <- c("Group", group_high, group_low)
  
  # --- 2.3 第一轮 DESeq2 (Raw Analysis) ---
  message("运行原始 DESeq2...")
  dds <- DESeqDataSetFromMatrix(countData = df_final, colData = coldata_sub, design = ~ Group)
  dds <- DESeq(dds)
  res_raw <- results(dds, contrast = contrast_vec) 
  
  # 保存差异分析结果
  res_raw_df <- as.data.frame(res_raw)
  res_raw_df$gene <- rownames(res_raw_df)
  write.table(res_raw_df, paste0(tissue, '_BW_DESeq2_raw.txt'), sep = '\t', row.names = F)
  
  # --- 2.4 SVA 批次效应分析 ---
  message("正在进行 SVA 分析...")
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  dat_for_sva <- assay(vsd)
  
  # [新增功能] 保存原始数据的 PCA 坐标表格
  # --------------------------------------------------------------------
  pca_raw <- prcomp(t(dat_for_sva))
  pca_raw_df <- as.data.frame(pca_raw$x) # 提取 PC1, PC2...
  # 合并 metadata 方便查看
  pca_raw_out <- cbind(Sample = rownames(pca_raw_df), 
                       coldata_sub[rownames(pca_raw_df), c("Group", "Project")], 
                       pca_raw_df)
  # 保存为表格
  write.table(pca_raw_out, paste0(tissue, "_PCA_Coordinates_Raw.txt"), sep = '\t', row.names = F, quote = F)
  message(paste("已保存原始 PCA 坐标表:", paste0(tissue, "_PCA_Coordinates_Raw.txt")))
  # --------------------------------------------------------------------

  mod <- model.matrix(~ Group, data = coldata_sub)
  mod0 <- model.matrix(~ 1, data = coldata_sub)
  
  n.sv <- num.sv(dat_for_sva, mod, method = "be")
  message(paste("检测到的 SV 数量:", n.sv))
  
  p1 <- plot_pca(dat_for_sva, coldata_sub, paste0(tissue, " - Raw")) + theme(legend.position = "none")
  p2 <- NULL 
  
  if (n.sv > 0) {
    svobj <- sva(dat_for_sva, mod, mod0, n.sv = n.sv)
    sv_data <- as.data.frame(svobj$sv)
    colnames(sv_data) <- paste0("SV", 1:ncol(sv_data))
    
    # 检查相关性
    r_sq_values <- sapply(1:n.sv, function(i) {
      fit <- lm(sv_data[, i] ~ as.numeric(coldata_sub$Group))
      summary(fit)$r.squared
    })
    
    bad_svs_idx <- which(r_sq_values > 0.5)
    
    if(length(bad_svs_idx) > 0) {
      message(paste("警告：剔除与分组高度相关的 SV:", paste(names(r_sq_values)[bad_svs_idx], collapse=",")))
      final_svs_names <- colnames(sv_data)[-bad_svs_idx]
    } else {
      final_svs_names <- colnames(sv_data)
    }
    
    if(length(final_svs_names) > 0) {
      coldata_with_sv <- cbind(coldata_sub, sv_data[, final_svs_names, drop=FALSE])
      design_formula <- as.formula(paste0("~ ", paste(c(final_svs_names, "Group"), collapse = " + ")))
      
      message("运行去批次后的 DESeq2...")
      dds_sva <- DESeqDataSetFromMatrix(countData = df_final, 
                                        colData = coldata_with_sv, 
                                        design = design_formula)
      dds_sva <- DESeq(dds_sva)
      
      res_sva <- results(dds_sva, contrast = contrast_vec)
      res_sva_df <- as.data.frame(res_sva)
      res_sva_df$gene <- rownames(res_sva_df)
      write.table(res_sva_df, paste0(tissue, '_BW_DESeq2_rmbatcheffect.txt'), sep = '\t', row.names = F)
      
      # 准备去批次后的数据
      sv_matrix_clean <- as.matrix(sv_data[, final_svs_names, drop=FALSE])
      dat_cleaned <- removeBatchEffect(dat_for_sva, 
                                       covariates = sv_matrix_clean, 
                                       design = mod)
      
      # [新增功能] 保存去批次后的 PCA 坐标表格
      # --------------------------------------------------------------------
      pca_clean <- prcomp(t(dat_cleaned))
      pca_clean_df <- as.data.frame(pca_clean$x)
      pca_clean_out <- cbind(Sample = rownames(pca_clean_df), 
                             coldata_sub[rownames(pca_clean_df), c("Group", "Project")], 
                             pca_clean_df)
      write.table(pca_clean_out, paste0(tissue, "_PCA_Coordinates_SVA_Corrected.txt"), sep = '\t', row.names = F, quote = F)
      message(paste("已保存去批次后 PCA 坐标表:", paste0(tissue, "_PCA_Coordinates_SVA_Corrected.txt")))
      # --------------------------------------------------------------------
      
      p2 <- plot_pca(dat_cleaned, coldata_sub, paste0(tissue, " - After SVA"))
      
    } else {
      message("剔除坏 SV 后无剩余 SV，跳过矫正步骤。")
      p2 <- plot_pca(dat_for_sva, coldata_sub, paste0(tissue, " - No Valid SVs"))
    }
    
  } else {
    message("没有检测到显著的批次效应，跳过矫正。")
    p2 <- plot_pca(dat_for_sva, coldata_sub, paste0(tissue, " - No SVs Found"))
  }
  
  # --- 2.5 保存图片 ---
  pdf(paste0('20260119_SVA_', tissue, '.pdf'), width = 12, height = 5)
  print(p1 + p2)
  dev.off()
  
  message(paste0(tissue, " 处理完成！\n"))
}

cat("所有器官处理完毕！")
