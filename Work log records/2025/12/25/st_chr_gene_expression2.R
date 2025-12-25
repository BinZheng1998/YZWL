# ==============================================================================
# 0. 环境设置与库加载
# ==============================================================================
library(Seurat)
library(dplyr)
library(tibble)
library(rtracklayer) # 用于读取GTF
library(ggplot2)
library(ggrepel)     # 用于标签防重叠
library(patchwork)   # 用于拼图
library(tidyr)

# 设置工作目录 (根据你的实际情况修改)
# setwd("~/project/analysis/")

# ==============================================================================
# 1. 数据准备函数：生成 Pseudo-bulk 矩阵
# ==============================================================================
# 说明：此函数将空间转录组数据按样本聚合成 Pseudo-bulk
get_pseudobulk <- function(rds_path, sample_name) {
  if (!file.exists(rds_path)) stop(paste("文件不存在:", rds_path))
  
  message(paste("正在处理:", sample_name, "from", rds_path))
  obj <- readRDS(rds_path)
  
  # 确保默认 Assay 正确
  DefaultAssay(obj) <- "Spatial"
  
  # 聚合 Count (Seurat V4/V5 兼容写法)
  pb_list <- AggregateExpression(obj, group.by = "orig.ident", assays = "Spatial", slot = "counts", return.seurat = FALSE)
  mat <- pb_list$Spatial
  
  # 重命名列名以防混淆
  colnames(mat) <- paste0(sample_name, "_", colnames(mat))
  
  # 转为数据框并保留基因名
  df <- as.data.frame(mat) %>% 
    rownames_to_column("gene_id")
  
  # 清理内存
  rm(obj, mat, pb_list)
  gc()
  return(df)
}

# ==============================================================================
# 2. 核心分析工具函数
# ==============================================================================

# --- 2.1 加载 GTF 并过滤常染色体 ---
load_sc_data <- function(cpm_matrix, gtf_path) {
  message(">> 正在读取 GTF 注释文件...")
  gtf <- import(gtf_path)
  
  # 整理基因坐标信息
  gene_map <- as.data.frame(gtf) %>%
    filter(type == "gene") %>%
    # 剔除 Z, W, MT 染色体 (只分析常染色体以避免剂量补偿干扰)
    filter(!seqnames %in% c("chrZ", "chrW", "chrMT", "Z", "W", "MT")) %>%
    dplyr::select(gene_id, chromosome = seqnames, any_of("gene_name")) %>% 
    distinct()
  
  message(">> 正在匹配矩阵与注释...")
  expr_mat <- as.matrix(cpm_matrix)
  
  # 取交集：矩阵中存在的基因 AND 常染色体上的基因
  common_genes <- intersect(rownames(expr_mat), gene_map$gene_id)
  
  # 过滤数据
  expr_mat_filtered <- expr_mat[common_genes, , drop = FALSE]
  gene_map_filtered <- gene_map %>% filter(gene_id %in% common_genes)
  
  message(sprintf("   原始基因数: %d", nrow(expr_mat)))
  message(sprintf("   分析用基因数 (常染色体): %d", length(common_genes)))
  
  return(list(expr_mat = expr_mat_filtered, gene_map = gene_map_filtered))
}

# --- 2.2 核心计算函数 (支持 Trim 和 No Trim) ---
analyze_chromosome_metrics <- function(expr_mat, gene_map, top_trim_prop = 0.05) {
  
  # 1. 预处理：Log2(CPM + 1)
  expr_mat_log <- log2(expr_mat + 1)
  
  # 2. 计算每个基因的均值
  gene_means <- rowMeans(expr_mat_log)
  
  # 3. 合并染色体信息
  gene_data <- data.frame(
    gene_id = rownames(expr_mat_log),
    expression = gene_means
  ) %>%
    inner_join(gene_map, by = "gene_id")
  
  # 4. 按染色体计算统计指标
  chrom_summary <- gene_data %>%
    group_by(chromosome) %>%
    summarise(
      n_total = n(),
      # 广度 (Breadth): 活跃基因比例 (CPM > 1)
      n_active = sum(expression > 1),
      ratio_active = n_active / n_total,
      
      # 强度 (Intensity): 计算截尾均值
      # 如果 top_trim_prop = 0，则 quantile 结果为最大值，即不修剪
      # 如果 top_trim_prop = 0.05，则去除前 5% 的极端高表达值
      avg_expression = mean(expression[expression <= quantile(expression, 1 - top_trim_prop)]),
      
      .groups = 'drop'
    )
  
  return(chrom_summary)
}

# --- 2.3 可视化函数 ---
plot_landscape <- function(chrom_summary, method_label) {
  
  # 1. 染色体分类与排序
  plot_data <- chrom_summary %>%
    mutate(chr_name = gsub("chr", "", chromosome)) %>%
    # 辅助排序列：将数字字符转为数字，非数字给极大值放在最后
    mutate(chr_num_sort = as.numeric(ifelse(grepl("^\\d+$", chr_name), chr_name, 999))) %>%
    arrange(chr_num_sort, chr_name) %>%
    mutate(chr_name = factor(chr_name, levels = unique(chr_name))) %>%
    # 定义类型
    mutate(type = case_when(
      chr_name %in% as.character(1:9) ~ "Macro (1-9)",
      chr_name %in% c(as.character(29:38), "16") ~ "Dot (29-38, 16)", 
      TRUE ~ "Micro (Others)"
    )) %>%
    mutate(type = factor(type, levels = c("Macro (1-9)", "Micro (Others)", "Dot (29-38, 16)")))
  
  # 2. 绘图
  p <- ggplot(plot_data, aes(x = ratio_active, y = avg_expression, color = type)) +
    geom_point(size = 5, alpha = 0.8) +
    geom_text_repel(aes(label = chr_name), size = 4, show.legend = FALSE, max.overlaps = 20) +
    
    # 辅助线
    geom_vline(xintercept = median(plot_data$ratio_active), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = median(plot_data$avg_expression), linetype = "dashed", color = "grey70") +
    
    scale_color_manual(values = c("#4dbbd5ff", "#e64b35ff", "#00a087ff")) +
    theme_classic(base_size = 14) +
    labs(
      x = "Breadth: Ratio of Active Genes (CPM > 1)",
      y = paste0("Intensity: ", ifelse(grepl("Trim", method_label), "Robust Mean", "Mean"), " log2(CPM+1)"),
      title = "Chromosome Transcriptional Landscape",
      subtitle = paste0("Method: ", method_label),
      color = "Chromosome Type"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
  
  return(p)
}

# ==============================================================================
# 3. 执行主流程 (Pipeline Execution)
# ==============================================================================

# --- 步骤 A: 生成并保存 CPM 矩阵 (如果已生成过 CSV，可跳过此步直接读取) ---
if (TRUE) { # 设置为 FALSE 可跳过此步
  message(">>> Step A: 生成 Pseudo-bulk 矩阵...")
  
  # 请确保路径正确
  df_E2.5 <- get_pseudobulk('../../00_data/E2.5_sp_new_seurat_v5_cluster.rds', "E2.5")
  df_E3.5 <- get_pseudobulk('../../00_data/E3.5_sp_new_seurat_v5_cluster.rds', "E3.5")
  df_E4.5 <- get_pseudobulk('../../00_data/E4.5_sp_new_seurat_v5_cluster.rds', "E4.5")
  
  # 合并
  final_df <- df_E2.5 %>%
    full_join(df_E3.5, by = "gene_id") %>%
    full_join(df_E4.5, by = "gene_id")
  final_df[is.na(final_df)] <- 0
  
  # 转换为矩阵
  final_mat <- final_df %>% column_to_rownames("gene_id") %>% as.matrix()
  
  # 保存原始 Count
  write.csv(final_mat, "chicken_embryo2pseudo_bulk_union_matrix.csv")
  
  # 计算 CPM 并保存
  # 增加安全性：防止列和为0导致NaN
  col_sums <- colSums(final_mat)
  col_sums[col_sums == 0] <- 1 
  cpm_mat <- t(t(final_mat) / col_sums * 1e6)
  
  # 清理中间变量
  rm(df_E2.5, df_E3.5, df_E4.5, final_df, final_mat)
  gc()
}

# --- 步骤 B: 加载数据与 GTF ---
# 定义 GTF 路径
gtf_file <- '~/project/00_ref/20250529_chicken.gtf' 

message(">>> Step B: 加载与过滤...")
sc_data <- load_sc_data(cpm_mat, gtf_file)

# --- 步骤 C: 运行两种分析模式 ---
message(">>> Step C: 运行分析...")

# 模式 1: No Trimming (保留所有数据)
message("   Running: No Trimming...")
res_no_trim <- analyze_chromosome_metrics(
  sc_data$expr_mat, 
  sc_data$gene_map, 
  top_trim_prop = 0 # 0% 修剪
)

# 模式 2: Robust / Trimmed (去除前 2.5% 高表达异常值)
message("   Running: Robust (Trim 5%)...")
res_trimmed <- analyze_chromosome_metrics(
  sc_data$expr_mat, 
  sc_data$gene_map, 
  top_trim_prop = 0.025 # 5% 修剪
)

# --- 步骤 D: 绘图与拼图 ---
message(">>> Step D: 生成图表...")

p1 <- plot_landscape(res_no_trim, "No Trimming (Raw Mean)")
p2 <- plot_landscape(res_trimmed, "Robust (Top 5% Outliers Removed)")

# 使用 patchwork 进行拼图
final_plot <- p1 + p2 + 
  plot_layout(guides = "collect") & # 共用图例
  theme(legend.position = "bottom")

# 打印查看
print(final_plot)

# 保存
ggsave("Comparison_Chromosome_Landscape_Trim_vs_NoTrim.pdf", final_plot, width = 16, height = 8)
ggsave("Robust_Only_Chromosome_Landscape.pdf", p2, width = 8, height = 8) # 单独保存 Robust 版

# 保存统计表格
write.csv(res_trimmed, "Chromosome_Metrics_Robust.csv", row.names = FALSE)
write.csv(res_no_trim, "Chromosome_Metrics_NoTrim.csv", row.names = FALSE)

message(">>> 全部完成！")
