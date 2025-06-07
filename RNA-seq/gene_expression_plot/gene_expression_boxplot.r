library(readxl)
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggsignif)

df <- read.csv('merged_gene_abund_results.csv', sep = ',', header = FALSE)
data <- read_excel('E:/07-转录组相关/01-鸡/rna/metadata/rnaseq_metadata.xlsx',sheet = 1)
data2 <- data[, c(2, 6, 9)]
colnames(data2) <- c('Sample', 'Breed', 'Tissue')

gene_list <- c('PLAG1', 'KHDRBS3', 'LCORL')  # 示例基因列表，可根据需要修改

results_df <- data.frame(Gene = character(), 
                         Mean_High_BW = numeric(), 
                         Mean_Low_BW = numeric(), 
                         SD_High_BW = numeric(), 
                         SD_Low_BW = numeric(), 
                         P_Value = numeric(), 
                         stringsAsFactors = FALSE)

for (gene in gene_list) {
  df2 <- df[df$V1 == paste0('gene-', gene) | df$V1 == 'Col1', ]
  df3 <- df2[, -c(2, 3)]
  df4 <- t(df3)
  colnames(df4) <- c('Sample', gene)
  
  final <- merge(df4, data2, by = 'Sample')
  final2 <- final[-grep(c('Adrenal gland|Beak|Jejunum|Colon|Pancreas|Shank|proventriculus|Egg|Gallbladder|Cochlea|Pulmonary Artery|Missing'), final$Tissue), ]
  final3 <- final2[final2[[gene]] > 0.00000001, ]
  final3[[gene]] <- log10(as.numeric(final3[[gene]]))
  
  final$BW <- NA
  final3$BW[final3$Breed %in% c("Baier Huang", "Huang shan Black", "Tibetan Chicken", "Red Jungle Fowl", "Fayoumi", "Chahua Chicken")] <- "Low BW"
  final3$BW[final3$Breed %in% c("White Plymouth Rock", "Rhode Island White", "Rhode Island Red", "Cornish", "Cobb", "Minorca", "Ross", "Ross 308", "Arbor Acres")] <- "High BW"
  
  final3_mean <- final3 %>% group_by(Tissue) %>% summarise(sd = sd(.data[[gene]]), value = mean(.data[[gene]]))
  final3_mean <- final3_mean[order(-final3_mean$value), ]
  final3_mean$Tissue <- factor(final3_mean$Tissue, levels = final3_mean$Tissue)
  
  final4 <- final3[final3$BW %in% c("Low BW", "High BW"), ]
  final4_mean <- final4 %>% group_by(BW) %>% summarise(sd = sd(.data[[gene]]), value = mean(.data[[gene]]))
  
  group1 <- final4[final4$BW == 'High BW', gene]
  group2 <- final4[final4$BW == 'Low BW', gene]
  test_result <- t.test(group1, group2)
  
  p_value <- test_result$p.value
  mean_high <- mean(group1, na.rm = TRUE)
  mean_low <- mean(group2, na.rm = TRUE)
  sd_high <- sd(group1, na.rm = TRUE)
  sd_low <- sd(group2, na.rm = TRUE)
  
  results_df <- rbind(results_df, data.frame(Gene = gene, 
                                             Mean_High_BW = mean_high, 
                                             Mean_Low_BW = mean_low, 
                                             SD_High_BW = sd_high, 
                                             SD_Low_BW = sd_low, 
                                             P_Value = p_value))
  
  p2 <- ggplot(final4_mean, aes(BW, value, fill = BW)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = 0.5) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 3)) +
    scale_fill_manual(values = c("black", "grey90")) +
    theme_bw() +
    ylab(gene) +  
    xlab('') +
    theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    geom_signif(data = final4, mapping = aes(BW, .data[[gene]]), comparisons = list(c("High BW", "Low BW")), map_signif_level = FALSE, tip_length = 0, test = 't.test')
  
  ggsave(p2, file = paste0(gene, '_High_vs_Low_BW.png'), dpi = 300, width = 2.5, height = 3)
}

write.table(results_df, file = "t_test_results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
