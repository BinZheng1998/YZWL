setwd('/Users/zhengbin/Desktop/崖州湾/群体进化/GO/Cattle/selection/top1/GO/')
library(dplyr)
library(clusterProfiler)
library(DOSE)
library(GO.db)
library(org.Bt.eg.db)
gene_files <- list.files(path = '../gene', pattern = "*genes.txt", full.names = TRUE)

# 循环处理每个基因文件
for (file in gene_files) {
  # 读取基因文件
  gene <- read.table(file, header = TRUE, sep = "\t")
  
  # 提取基因列表
  genes_list <- strsplit(gene[[5]], ",")
  all_genes <- unlist(genes_list)
  unique_genes <- unique(all_genes)
  symbol <- as.character(unique_genes)
  
  # 转换基因符号为 ENTREZID
  eg <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Bt.eg.db")
  id <- as.character(eg[, 2])
  
  # 进行 GO 富集分析
  ego <- enrichGO(gene = id,
                  OrgDb = org.Bt.eg.db,
                  ont = "all",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1,
                  qvalueCutoff = 1,
                  readable = TRUE)
  
  # 过滤 p 值小于 0.05 的结果
  ego <- as.data.frame(ego)
  ego <- ego[ego$pvalue < 0.05, ]
  
  # 生成输出文件名
  output_file <- sub("genes.txt", "GO.txt", basename(file))
  
  # 保存结果
  write.table(ego, file = output_file, sep = "\t", row.names = FALSE)
}
