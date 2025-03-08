library(Seurat)
library(pheatmap)
df <- readRDS("input.rds")
Idents(df)<- df$subcluster
#如果是单细胞数据则是av.exp<- AverageExpression(df)$RNA
av.exp<- AverageExpression(df)$Spatial
features=names(tail(sort(apply(av.exp, 1, sd)),2000))
av.exp <- as.matrix(av.exp)
av.exp <- cor(av.exp, method= "spearman")
p<-pheatmap::pheatmap(av.exp)
ggsave(p,file="E4.5.pdf")
