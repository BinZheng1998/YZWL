setwd("/home/bzheng/project/04_rnaseq/05_plot/")
# 加载所需的包
library(readxl)
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggsignif)

# 读取数据（只需执行一次，放在循环外）
df <- read.csv('../04_expression_result/merged_gene_abund_results.csv', sep = ',', header = FALSE)
data <- read.table("../metadata/chicken_rnaseq_metadata.txt",sep = '\t',header = T,fill=T)
data2 <- data[, c(2, 9, 12)]
colnames(data2) <- c('Sample', 'BW', 'Tissue')
data2 <- data2[data2$Tissue == "Liver",]
gene <- read.table("~/project/01_evolution/convergent_evolution/07_result/20250620_homo_gene/BW/5species_pos_BW_top5_50kb_homo_genes.txt",sep = '\t',header = T)
gene_list <- gene$chicken_id
results_df <- data.frame(Gene = character(), 
                         Mean_High_BW = numeric(), 
                         Mean_Low_BW = numeric(), 
                         SD_High_BW = numeric(), 
                         SD_Low_BW = numeric(), 
                         DC = numeric(),
                         P_Value = numeric(), 
                         stringsAsFactors = FALSE)

# 开始循环分析
for (gene in gene_list) {
  # 选择当前基因的数据
  df2 <- df[df$V1 == gene | df$V1 == 'Col1', ]
  df3 <- df2[, -c(2, 3)]
  df4 <- t(df3)
  colnames(df4) <- c('Sample', gene)
  
  # 合并数据
  final <- merge(df4, data2, by = 'Sample')
  final <- final[!grepl("Middle|N/A", final$BW), ]
  final2 <- final[final[[gene]] > 0.00000001, ]
  if (nrow(final2) < 6) {
    message(sprintf("⚠️ 基因 %s 过滤后样本少于6个", gene))
    next
  }
  # 如果 High 或 Low 组样本数量 < 3 则跳过
  group_counts <- table(final2$BW[final2$BW %in% c("Low", "High")])
  low_count <- as.numeric(group_counts["Low"] %||% 0)
  high_count <- as.numeric(group_counts["High"] %||% 0)
  if (low_count < 3 || high_count < 3) {
    message(sprintf("⚠️ 基因 %s 在 High (%d) 或 Low (%d) 中样本数少于3，跳过分析", gene, high_count, low_count))
    next
  }
  
  final2[[gene]] <- log10(as.numeric(final2[[gene]]))
  
  # 计算每个组织的均值和标准差
  final2_mean <- final2 %>% group_by(Tissue) %>% summarise(sd = sd(.data[[gene]]), value = mean(.data[[gene]]))
  final2_mean <- final2_mean[order(-final2_mean$value), ]
  final2_mean$Tissue <- factor(final2_mean$Tissue, levels = final2_mean$Tissue)
  
  # 体重分析
  final3 <- final2[final2$BW %in% c("Low", "High"), ]
  final3_mean <- final3 %>% group_by(BW) %>% summarise(sd = sd(.data[[gene]]), value = mean(.data[[gene]]))
  
  # 进行 t 检验并提取结果
  group1 <- final3[final3$BW == 'High', gene]
  group2 <- final3[final3$BW == 'Low', gene]
  test_result <- t.test(group1, group2)
  
  # 提取统计信息
  p_value <- test_result$p.value
  mean_high <- mean(group1, na.rm = TRUE)
  mean_low <- mean(group2, na.rm = TRUE)
  sd_high <- sd(group1, na.rm = TRUE)
  sd_low <- sd(group2, na.rm = TRUE)
  FC_highvslow <- mean_high/mean_low
  
  # 将结果添加到数据框中
  results_df <- rbind(results_df, data.frame(Gene = gene, 
                                             Mean_High_BW = mean_high, 
                                             Mean_Low_BW = mean_low, 
                                             SD_High_BW = sd_high, 
                                             SD_Low_BW = sd_low,
                                             FC = FC_highvslow,
                                             P_Value = p_value))
  
  y_max <- max(final3[[gene]], na.rm = TRUE) * 1.3
  
  # 绘制图表
  p2 <- ggplot(final3, aes(BW, .data[[gene]], fill = BW)) +
    geom_violin(color="white") + 
    geom_boxplot(width=0.1)+
    scale_fill_manual(values = c("#B2CDE2", "#4A81B1")) +
    theme_bw() +
    ylab(gene) +  # 使用当前基因名作为y轴标签
    xlab('') +
    scale_y_continuous(limits = c(NA, y_max)) +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    geom_signif(data = final3, mapping = aes(BW, .data[[gene]]), comparisons = list(c("High", "Low")), 
                vjust = -0.2,y_position = y_max * 0.9,
                map_signif_level = F, tip_length = 0, test = 't.test')
  
  # 保存图表，文件名包含当前基因名
  ggsave(p2, file = paste0(gene, '_High_vs_Low_Liver_BWpos.pdf'), dpi = 300, width = 2.5, height = 3)
}
results_df$P.adjust <- p.adjust(results_df$P_Value)
results_df <- unique(results_df)
# 在循环结束后，将数据框写入 txt 文件
write.table(results_df, file = "BW_pos_top5_50kb_gene_Liver_ttest_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)


