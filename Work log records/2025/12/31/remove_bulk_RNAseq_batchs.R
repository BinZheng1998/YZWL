setwd('~/project/04_rnaseq/chicken/05_result/251230_PCA/')
library(DESeq2)
library(ggplot2)
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
df <- df[df$Filter == "PASS",]
df2 <- df[, c(2, 4,8,11, 14,16)]
head(df2)
colnames(df2) <- c('Sample','Project','Breed', 'BW', 'Tissue','Instrument')
df2 <- df2[df2$Tissue %in% c("Liver","Muscle","Hypothalamus","Adipose"),]
df2 <- df2[df2$BW %in% c("High","Low"),]
df2$Group <- paste0(df2$BW,sep = '_',df2$Tissue)
df3 <- df2[df2$Tissue == "Hypothalamus",]

valid_cols <- intersect(df3$Sample, colnames(data4))
data5 <- data4[, valid_cols]
keep <- rowSums(data5 > 5) >= 10 
data5 <- data5[keep, ]
data6 <- data5 +1 

head(data6)
coldata <- df3[df3$Sample %in% colnames(data6),]
rownames(coldata) <- coldata[,1]
coldata$Group <- as.factor(coldata$Group)
dds<-DESeqDataSetFromMatrix(countData=data6,colData=coldata,design=~Group)
vsd <- vst(dds, blind=FALSE)
#dds<-DESeq(dds)
pcaData <- plotPCA(vsd, intgroup=c("Group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=Group)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  ggtitle("PCA Plot of RNA-seq Samples")
 
library(sva)
coldata$Project <- as.factor(coldata$Project)
coldata$Breed <- as.factor(coldata$Breed)
coldata$Instrument <- as.factor(coldata$Instrument)
coldata$Group <- as.factor(coldata$Group)

library(limma)

dds_simple <- DESeqDataSetFromMatrix(countData = data6, colData = coldata, design = ~ 1)
vsd_simple <- vst(dds_simple, blind = FALSE)
mat <- assay(vsd_simple)

design_mat <- model.matrix(~ Group, data = coldata)
adj_mat <- removeBatchEffect(mat, 
                             batch = coldata$Breed, 
                             batch2 = coldata$Instrument, 
                             design = design_mat)

pca_res <- prcomp(t(adj_mat))
pca_df <- as.data.frame(pca_res$x)
percentVar_adj <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))

pca_df$Group <- coldata$Group
pca_df$Breed <- coldata$Breed
ggplot(pca_df, aes(PC1, PC2, color = Breed, shape = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = rainbow(13)) + 
  scale_shape_manual(values = c(16, 17)) + 
  xlab(paste0("PC1: ", percentVar_adj[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_adj[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA after Forced Project Correction (limma)")



# 尝试保留最重要的两个协变量：Project 和 Breed
# 如果依然报错，请继续删减 Breed，只保留 ~ Project + Group
dds_final <- DESeqDataSetFromMatrix(countData = data6, 
                                    colData = coldata, 
                                    design = ~ Project + Group)

# 运行差异分析
dds_final <- DESeq(dds_final)

# 提取结果
res <- results(dds_final, contrast = c("Group", "High_Hypothalamus", "Low_Hypothalamus"))
res <- as.data.frame(res)
summary(res)
res1 <- res
res1$gene <- rownames(res1)
res2 <-  res1[!grepl("STRG.", res1$gene), ]
