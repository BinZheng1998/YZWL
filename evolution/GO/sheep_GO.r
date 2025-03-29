setwd('E:/02-群体进化/多物种/GO/202503/sheep/selection_BW/top1/GO/')
library(dplyr)
library(clusterProfiler)
library(DOSE)
library(GO.db)
require(AnnotationHub)
hub <- AnnotationHub::AnnotationHub()
query(hub,"Ovis aries")
sheep <- hub[['AH114632']]
length(keys(sheep))
columns(sheep)


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
  symbol <- gsub(" ","",symbol)
  
  # 转换基因符号为 ENTREZID
  eg <- bitr(symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = sheep)
  id <- as.character(eg[, 2])
  
  # 进行 GO 富集分析
  ego <- enrichGO(gene = id,
                  OrgDb = sheep,
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
