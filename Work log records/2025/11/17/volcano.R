setwd('~/project/04_rnaseq/chicken/05_result/251113_volcano/')
library(DESeq2)
data <- read.csv('../../03_stringtie_2.2.1/re_quant/gene_count_matrix.csv')
data$gene_id <- gsub(".*\\|", "", data$gene_id)
data <- as.data.frame(data)
data1 <- data[!duplicated(data$gene_id), ]
rownames(data1) <- data1$gene_id
data2 <- data1[,-1]
data3 <- data2[, colSums(data2, na.rm = TRUE) >= 1000000]

detected_genes_gt1 <- colSums(data3 > 1)
summary(detected_genes_gt1)
hist(detected_genes_gt1, 
     main = "Histogram of Detected Genes per Sample",
     xlab = "Number of Detected Genes (TPM > 1)",
     breaks = 30) 
#去除极端值
Q1 <- quantile(detected_genes_gt1, 0.25)
Q3 <- quantile(detected_genes_gt1, 0.75)
IQR_val <- IQR(detected_genes_gt1)
whisker_lower <- Q1 - 1.5 * IQR_val
whisker_upper <- Q3 + 1.5 * IQR_val
samples_to_keep <- names(detected_genes_gt1[detected_genes_gt1 >= whisker_lower & detected_genes_gt1 <= whisker_upper])
data4 <- data3[, samples_to_keep]


df <- read.table("../../metadata/chicken_rnaseq_metadata.txt",sep = '\t',header = T,fill=T)
df2 <- df[, c(2, 11, 14)]
colnames(df2) <- c('Sample', 'BW', 'Tissue')
df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose"),]
df2 <- df2[df2$BW %in% c("High","Low"),]
df2$Group <- paste0(df2$BW,sep = '_',df2$Tissue)
df3 <- df2[df2$Tissue == "Hypothalamus",]

valid_cols <- intersect(df3$Sample, colnames(data4))
data5 <- data4[, valid_cols]
data6 <- data5 +1 

head(df3)
coldata <- df3[df3$Sample %in% colnames(data6),]
rownames(coldata) <- coldata[,1]
coldata$Group <- as.factor(coldata$Group)
dds<-DESeqDataSetFromMatrix(countData=data6,colData=coldata,design=~Group)
dds<-DESeq(dds)
res<-results(dds)
res<-as.data.frame(res)

genes1 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS/BW_convergent/4species_intersect_pos_BW_top5_50kb_homo_genes.txt',header = T)
genes1$type <- 'pos'
genes2 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS/BW_convergent/4species_intersect_neg_BW_top5_50kb_homo_genes.txt',header = T)
genes2$type <- 'neg'
genes3 <- rbind(genes1,genes2)
genes <- intersect(genes1$chicken, genes2$chicken)
genes

genes4 <- genes3[,c(5,6)]
genes4$type[genes4$chicken %in% genes] <- 'both'
genes4 <- unique(genes4)

res1 <- res[rownames(res) %in% genes4$chicken,]

match_indices <- match(rownames(res1), genes4$chicken)
res1$Type <- genes4$type[match_indices]

library(ggrepel)
library(ggplot2)
head(res1)
pthreshold <- 0.05
fcthreshold <- 0.5 

sig_genes_df <- res1[res1$padj < pthreshold & abs(res1$log2FoldChange) > fcthreshold, ]
if (nrow(sig_genes_df) > 0) {
  top5_genes <- rownames(sig_genes_df[order(sig_genes_df$padj), ])[1:min(5, nrow(sig_genes_df))]
  res1$label <- ifelse(rownames(res1) %in% top5_genes, as.character(rownames(res1)), "")
} else {
  res1$label <- ""
}


up_count <- sum(res1$log2FoldChange > fcthreshold & res1$padj < pthreshold)
down_count <- sum(res1$log2FoldChange < -fcthreshold & res1$padj < pthreshold)
stable_count <- nrow(res1) - up_count - down_count

pdf('chicken_Hypothalamus_volcano_4species_top5_neg_pos.pdf',width = 8,height = 6)
ggplot(data = res1,aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = log2FoldChange, size = -log10(padj),shape=Type)) + 
  #scale_shape_manual(values = c("pos" = 17, "neg" = 16, "both" = 18)) +
  scale_color_gradient2(low = "blue", 
                        mid = "grey", 
                        high = "red",
                        midpoint = 0) + 
  scale_size_continuous(range = c(1, 8)) +
  scale_x_continuous(limits = c(-7,7))+
  scale_y_continuous(limits = c(0,250))+
  geom_text_repel(aes(label = label), 
                  size = 5,
                  box.padding = unit(0.6, "lines"),
                  point.padding = unit(0.7, "lines"),
                  segment.color = "black",
                  show.legend = FALSE,
                  max.overlaps = 100) +
  geom_vline(xintercept = c(-fcthreshold, fcthreshold), lty = 2, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(pthreshold), lty = 2, col = "black", lwd = 0.8) +
  theme_bw() +
  labs(x = "log2(FoldChange)", y = "-log10(Padj)", 
       #title = "Volcano Plot of Different Expression Proteins",
       color = "log2(FoldChange)", 
       size = "-log10(Padj)") + 
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "right")+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  annotate("segment",
           x = fcthreshold, y = 240, 
           xend = 6, yend = 240,
           arrow = arrow(length = unit(0.3, "cm")), 
           color = "red", lwd = 1) +
  annotate("segment",
           x = -fcthreshold, y = 240, 
           xend = -6, yend = 240,
           arrow = arrow(length = unit(0.3, "cm")),
           color = "blue", lwd = 1) +
  annotate("text",
           x = 4, 
           y = 250, 
          color = 'red',
           label = paste("Up:", up_count),,
           size = 5, fontface = "bold") +
  annotate("text",
           x = -4, 
           y = 250,
           color = 'blue',
           label = paste("Down:", down_count),
           size = 5, fontface = "bold") +
  annotate("text",
           x = 0, 
           y = 250,
           color = 'black',
           label = stable_count,
           size = 5, fontface = "bold") +
  annotate("text",
           x = 5, 
           y = 10,
           color = 'black',
           label = "bolditalic(P) == bold('0.05')",
           parse = TRUE,
           size = 5, fontface = "bold") 
dev.off()
