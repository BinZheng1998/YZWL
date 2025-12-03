setwd('~/project/04_rnaseq/sheep/2.result/251202_sheep_volcano/')
library(ggplot2)
library(DESeq2)
data <- read.table('../../1.data/SheepGTEx_v0_2512.Gene.count.txt',sep = '\t',header = T)
rownames(data) <- data$Geneid
colnames(data)
sample <- data.frame(colnames(data))

library(digest)
fingerprints <- apply(data, 2, digest)
dup_status <- duplicated(fingerprints)
expr_clean <- data[, !dup_status]
#ample <- data.frame(colnames(expr_clean))

metadata <- read.table('../../1.data/metadata/sheepGTEx_metadata.txt',sep = '\t',header = T,quote = "",fill = TRUE,comment.char = "")
metadata <- metadata[metadata$Filter == 'pass',]
metadata1 <- metadata[metadata$Tissue.cell.type == "Muscle",]
metadata1 <- metadata[metadata$Tissue.cell.type == "Liver",]
metadata1 <- metadata[metadata$Tissue.cell.subtype == "Hypothalamus",]
metadata1 <- metadata[metadata$Tissue.cell.type == "Adipose",]
#table(metadata1$Tissue.cell.subtype)
metadata1 <- metadata[metadata$Tissue.cell.subtype == c("Intermuscular adipose","Subcutaneous adipose","Perirenal adipose","Adipose"),]
metadata2 <- metadata1[,c("Sample","Body.weight.Group")]
metadata3 <- metadata2[metadata2$Body.weight.Group %in% c("low","high"),]


valid_cols <- intersect(metadata3$Sample, colnames(expr_clean))
#res <- setdiff(metadata3$sample1, colnames(expr_clean))
expr_clean1 <- expr_clean[, valid_cols]
expr_clean2 <- expr_clean1 +1 

coldata <- metadata3[metadata3$Sample %in% colnames(expr_clean2),]
rownames(coldata) <- coldata[,1]
coldata$Body.weight.Group <- as.factor(coldata$Body.weight.Group)
dds<-DESeqDataSetFromMatrix(countData=expr_clean2,colData=coldata,design=~Body.weight.Group)
dds<-DESeq(dds)
#这里将low BW 作为参考组，high BW作为实验组
res<-results(dds,contrast = c("Body.weight.Group","high","low"))#  参数的格式是：c("分组变量名", "测试组", "参考组")
res<-as.data.frame(res)
res$gene <- rownames(res)
write.table(res,'sheep_Hypothalamus_DESeq2.txt',sep = '\t',row.names = F)

genes1 <- read.table('~/project/01_evolution/convergent_evolution/07_result/20251123_BW_new_analysis/tables/4species_BW_convergent_genes.txt',sep = '\t',header = T)
match_indices <- match(rownames(res), genes1$sheep)
res$Type <- genes1$Selection[match_indices]
res1 <- na.omit(res)


library(ggrepel)
library(ggplot2)
head(res1)
pthreshold <- 0.05
fcthreshold <- 0.5

up_genes <- subset(res1, padj < pthreshold & log2FoldChange > fcthreshold)
up_genes <- up_genes[order(up_genes$padj), ]
top5_up <- head(rownames(up_genes), 5)
down_genes <- subset(res1, padj < pthreshold & log2FoldChange < -fcthreshold)
down_genes <- down_genes[order(down_genes$padj), ]
top5_down <- head(rownames(down_genes), 5)
top_genes_to_label <- c(top5_up, top5_down)
res1$label <- ifelse(rownames(res1) %in% top_genes_to_label,rownames(res1), "")

res1 <- na.omit(res1)
up_count <- sum(res1$log2FoldChange > fcthreshold & res1$padj < pthreshold)
down_count <- sum(res1$log2FoldChange < -fcthreshold & res1$padj < pthreshold)
stable_count <- nrow(res1) - up_count - down_count

X_MIN <- -7
X_MAX <- 7
Y_MIN <- 0
Y_MAX <- 15
X_RANGE <- X_MAX - X_MIN
Y_RANGE <- Y_MAX - Y_MIN

pdf('sheep_Liver_volcano_4species_top5_neg_pos.pdf',width = 8,height = 6)
ggplot(data = res1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = log2FoldChange, size = -log10(padj), shape = Type)) + 
  #scale_shape_manual(values = c("pos" = 17, "neg" = 16, "both" = 18)) +
  scale_color_gradient2(low = "blue", 
                        mid = "grey", 
                        high = "red",
                        midpoint = 0) + 
  scale_size_continuous(range = c(1, 8)) +
  scale_x_continuous(limits = c(X_MIN, X_MAX)) +
  scale_y_continuous(limits = c(Y_MIN, Y_MAX)) +
  geom_text_repel(aes(label = label), 
                  size = 5,fontface = "italic",
                  box.padding = unit(0.6, "lines"),
                  point.padding = unit(0.7, "lines"),
                  segment.color = "black",
                  show.legend = FALSE,
                  max.overlaps = 100) +
  geom_vline(xintercept = c(-fcthreshold, fcthreshold), lty = 2, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(pthreshold), lty = 2, col = "black", lwd = 0.8) +
  theme_bw() +
  labs(
  x = expression(log[2]("FoldChange")),
  y = expression(-log[10]("Padj")), 
       color = expression(log[2]("FoldChange")), 
       size = expression(-log[10]("Padj")))+ 
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 13),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        legend.text = element_text(size = 11), legend.title = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        legend.position = "right") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  annotate("segment",
           x = fcthreshold, 
           y = Y_MIN + (Y_RANGE * 0.92), 
           xend = X_MIN + (X_RANGE * 0.85), 
           yend = Y_MIN + (Y_RANGE * 0.92),
           arrow = arrow(length = unit(0.3, "cm"),type = "closed"), 
           color = "red", lwd = 1) +
  annotate("segment",
           x = -fcthreshold, 
           y = Y_MIN + (Y_RANGE * 0.92), 
           xend = X_MIN + (X_RANGE * 0.15), 
           yend = Y_MIN + (Y_RANGE * 0.92),
           arrow = arrow(length = unit(0.3, "cm"),type = "closed"),
           color = "blue", lwd = 1) +
  annotate("text",
           x = X_MIN + (X_RANGE * 0.7), 
           y = Y_MIN + (Y_RANGE * 0.96), 
           color = 'red',
           label = paste("BW High:", up_count),
           size = 5, fontface = "bold") +
  annotate("text",
           x = X_MIN + (X_RANGE * 0.3), 
           y = Y_MIN + (Y_RANGE * 0.96),
           color = 'blue',
           label = paste("BW Low:", down_count),
           size = 5, fontface = "bold") +
  annotate("text",
           x = X_MIN + (X_RANGE * 0.50), 
           y = Y_MIN + (Y_RANGE * 0.96),
           color = 'black',
           label = stable_count,
           size = 5, fontface = "bold") +
  annotate("text",
           x = X_MIN + (X_RANGE * 0.875), 
           y = 0,
           color = 'black',
           label = "bolditalic(P) == bold('0.05')",
           parse = TRUE,
           size = 5, fontface = "bold") 
dev.off()
