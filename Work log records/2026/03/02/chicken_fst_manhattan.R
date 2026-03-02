setwd('~/project/01_evolution/convergent_evolution/07_result/20260228_BW_new_result/plot/')

library(dplyr)
library(ggplot2)
library(readr)     # 读取速度快
library(patchwork) # 拼图专用

# ============================================================
# 1. 定义一个通用的绘图函数
#    input_file: 文件路径
#    y_label: Y轴标题
#    show_x_text: 是否显示X轴标签 (拼图时通常只在最底下一张图显示)
# ============================================================
plot_fst_window <- function(input_file, y_label, show_x_text = FALSE) {
  
  # 读取数据
  df <- read_tsv(input_file, show_col_types = FALSE)
  
  # 数据清洗
  # VCFtools输出通常列为: CHROM, BIN_START, BIN_END, N_VARIANTS, WEIGHTED_FST, MEAN_FST
  # 这里假设第6列是你需要的Fst值 (原代码中你取的是第6列)
  colnames(df)[6] <- "Fst_Value" 
  
  Fst_pos <- df %>%
    mutate(
      CHR = gsub("chr", "", CHROM),
      POS = as.numeric(BIN_START),
      Fst = as.numeric(Fst_Value)
    ) %>%
    filter(CHR %in% as.character(1:38)) %>%  # 筛选 1-38 号染色体
    mutate(CHR = as.numeric(CHR)) %>%         # 转为数值以正确排序
    arrange(CHR, POS)
  
  # 计算累积位置 (BPcum)
  don <- Fst_pos %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POS), .groups = 'drop') %>%
    mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    left_join(Fst_pos, ., by = "CHR") %>%
    arrange(CHR, POS) %>%
    mutate(
      BPcum = POS + total,
      color_group = ifelse(CHR %% 2 == 1, "A", "B") # 奇偶染色体不同色
    )
  
  # 计算 X 轴标签位置 (染色体中心)
  X_axis <- don %>%
    group_by(CHR) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = 'drop')
  
  # 处理标签：为了防拥挤，隐藏部分染色体编号
  hidden_chr <- c(11, 13, 15, 16, 17, 19, 20, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35, 37, 38)
  X_axis <- X_axis %>%
    mutate(label = ifelse(CHR %in% hidden_chr, "", as.character(CHR)))
  
  # 绘图
  my_colors <- c("A" = "grey50", "B" = "grey80")
  
  p <- ggplot(don, aes(x = BPcum, y = Fst, color = color_group)) +
    geom_point(size = 1, alpha = 0.8) +
    scale_color_manual(values = my_colors) +
    scale_x_continuous(label = X_axis$label, 
                       breaks = X_axis$center, 
                       expand = expansion(mult = c(0.01, 0.01))) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      axis.title = element_text(size = 18),
      axis.text.y = element_text(size = 15, color = "black"),
      legend.position = "none",
      axis.title.x = element_blank() # 默认移除X轴标题
    ) +
    ylab(y_label)
  
  # 如果不显示X轴文字 (针对 p1, p2, p3)
  if (!show_x_text) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    # 针对 p4 (最底部) 显示文字
    p <- p + theme(axis.text.x = element_text(size = 15, color = "black"))
  }
  
  return(p)
}

# ============================================================
# 2. 调用函数生成 4 张图
# ============================================================

# 定义文件路径
f1 <- '~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/chicken_5kwindow_2kstep_Fst.windowed.weir.fst'
f2 <- '~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/chicken_10kwindow_5kstep_Fst.windowed.weir.fst'
f3 <- '~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/chicken_20kwindow_10kstep_Fst.windowed.weir.fst'
f4 <- '~/project/01_evolution/convergent_evolution/01_chicken/new-res/selection_analysis_BW/DP4_2allele_50miss_maf005/chicken_50kwindow_25kstep_Fst.windowed.weir.fst'

# 生成图片
# p1, p2, p3 不显示 x轴标签，只有 p4 显示
cat("正在绘制 p1 (5k)...\n")
p1 <- plot_fst_window(f1, "Fst 5kb window", show_x_text = FALSE)

cat("正在绘制 p2 (10k)...\n")
p2 <- plot_fst_window(f2, "Fst 10kb window", show_x_text = FALSE)

cat("正在绘制 p3 (20k)...\n")
p3 <- plot_fst_window(f3, "Fst 20kb window", show_x_text = FALSE)

cat("正在绘制 p4 (50k)...\n")
p4 <- plot_fst_window(f4, "Fst 50kb window", show_x_text = TRUE)


# ============================================================
# 3. 拼图并保存
# ============================================================
cat("正在保存 PDF...\n")

pdf('chicken_fst.pdf', width = 24, height = 16)
# 使用 patchwork 语法拼图
p1 / p2 / p3 / p4
dev.off()

cat("完成！\n")
