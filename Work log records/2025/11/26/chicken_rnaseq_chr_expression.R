setwd('~/project/04_rnaseq/chicken/05_result/251126_chr_gene_expression/')
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)
library(rtracklayer)

exprt <- read.csv('../../04_expression/Matrix_TPM.csv')
exprt2 <- exprt[!duplicated(exprt$Geneid),]
rownames(exprt2) <- exprt2$Geneid
colnames(exprt2)[1] <- 'gene_id'

#gtf <- read.table('../../01_ref/rnaseq_ref/20250529_chicken.gtf',sep = '\t')
gtf <- import('../../01_ref/rnaseq_ref/20250529_chicken.gtf')
gtf_df <- as.data.frame(gtf)
gene_map <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(gene_id, chromosome = seqnames) %>% 
  distinct()
#write.table(gene_map[1],'chicken_gene_list.txt',row.names = F,col.names = F,quote = F)

merged_data <- inner_join(exprt2, gene_map, by = "gene_id")
chrom_tpm_result <- merged_data %>%
  group_by(chromosome) %>%
  # 关键点：这里把 sum 改成 mean
  summarise(across(where(is.numeric), mean)) %>% 
  column_to_rownames("chromosome")

final_result <- chrom_tpm_result[rownames(chrom_tpm_result) %in% paste0('chr',c(1:38)),]
sample_totals <- colSums(final_result)
keep_mask <- sample_totals >= 100
removed_samples <- names(sample_totals)[!keep_mask]
if(length(removed_samples) > 0) {
  cat("已剔除以下 TPM 总和过低的样本：\n")
  print(removed_samples)
} else {
  cat("所有样本的 TPM 总和均大于 100，无需剔除。\n")
}
final_result_filtered <- final_result[, keep_mask]




plot_data_box <- plot_data_box %>%
  mutate(type = case_when(
    chromosome %in% c(paste0("chr", 1:9), "chrZ", "chrW") ~ "Macro",
    chromosome %in% c(paste0("chr", 29:38),"chr16") ~ "Dot",
    TRUE ~ "Micro"
  )) %>%
  mutate(chromosome = gsub("chr", "", chromosome))

plot_data_box$chromosome <- reorder(plot_data_box$chromosome, plot_data_box$log_tpm, FUN = median)
plot_data_box$type <- factor(plot_data_box$type, levels = c("Macro", "Micro", "Dot"))

pdf("chicken_rnaseq_chr_expression_boxplot.pdf",height = 6,width = 12)
ggplot(plot_data_box, aes(x = chromosome, y = log_tpm)) +
  geom_boxplot(outlier.shape = NA, aes(color = type), fill = "white", width = 0.7) +
  # geom_jitter(width = 0.2, aes(color = type), size = 1.5, alpha = 0.7) +
  scale_color_manual(values = c(
    "Macro" = "#4E84C4", 
    "Micro" = "#F28E79", 
    "Dot"   = "#D14949"
  )) +
  
  theme_classic() +
  labs(color = "", y = "log2(TPM + 1)", x = "Chromosome") +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12)
  )
dev.off()
