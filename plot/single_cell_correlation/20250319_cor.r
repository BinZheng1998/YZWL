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
