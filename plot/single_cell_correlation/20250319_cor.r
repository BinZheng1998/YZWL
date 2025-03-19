library(Seurat)
df1 <- readRDS("df1.rds")
#df1 <- Read10X(data.dir="input1/",gene.column=1)
#df1 <- CreateSeuratObject(df1)

df2 <- readRDS("df2.rds")
#df2 <- Read10X(data.dir="input2/",gene.column=1)
#df2 <- CreateSeuratObject(df2)

df1$test <- "No_coverged_gene"
df2$test <- "Coverged_gene"

expr_matrix1 <- df1@assays$RNA@layer$data
expr_matrix2 <- df2@assays$RNA@layer$data

rownames(expr_matrix1) <- rownames(df1)
rownames(expr_matrix2) <- rownames(df2)
common_genes <- intersect(rownames(expr_matrix1), rownames(expr_matrix2))

avg_expr1 <- rowMeans(expr_matrix1[common_genes, ])
avg_expr2 <- rowMeans(expr_matrix2[common_genes, ])

cor_res <- cor.test(avg_expr1,avg_expr2,method="pearson")

cor_res$p.value
cor_res$estimate

avg_expr1 <- as.matrix(avg_expr1)
avg_expr2 <- as.matrix(avg_expr2)
df <- cbind(avg_expr1,avg_expr2)
colnames(df) <- c("avg_expr1","avg_expr2")
library(ggpubr)
p<-ggplot(df,aes(x=avg_expr1,y=avg_expr2))+
  geom_point(size = 3)+
  geom_smooth(method = 'lm',color='blue',formula = y ~ x)+
  stat_cor(method = 'spearman')+
  theme_bw()+
  xlab('Dataset1 Expression (log-normalized)')+
  ylab('Dataset2 Expression (log-normalized)')
ggsave(filename="cor.png",plot=p,dpi=500,height=5,width=5)
