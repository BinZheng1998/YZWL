setwd("~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS/BW_convergent/")
library(dplyr)
chicken <- read.table("../chicken/chicken_BW_neg_DCMS_top_001_merged_regions_50kb_genes.txt",sep = '\t',header = T)
genes_list <- strsplit(chicken[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
chicken_genes <- unique(symbol)

cattle <- read.table("../cattle/cattle_BW_neg_DCMS_top_001_merged_regions_50kb_genes.txt",sep = '\t',header = T)
genes_list <- strsplit(cattle[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
cattle_genes <- unique(symbol)

pig <- read.table("../pig/pig_BW_neg_DCMS_top_001_merged_regions_50kb_genes.txt",sep = '\t',header = T)
genes_list <- strsplit(pig[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
pig_genes <- unique(symbol)

sheep <- read.table("../sheep/sheep_BW_neg_DCMS_top_001_merged_regions_50kb_genes.txt",sep = '\t',header = T)
genes_list <- strsplit(sheep[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
sheep_genes <- unique(symbol)

dog <- read.table("../dog/dog_BW_neg_DCMS_top_001_merged_regions_50kb_genes.txt",sep = '\t',header = T)
genes_list <- strsplit(dog[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
dog_genes <- unique(symbol)


homogene <- read.table("~/project/10_RER/result/RER/RERconverge_res_20250708.txt",sep = '\t',header = T)
homogene <- homogene[,c(9,10,11,12,13)]
homogene[homogene == "-"] <- NA
homogene <- na.omit(homogene)
filtered_data <- homogene %>%
  filter(
    chicken %in% chicken_genes,
    cattle %in% cattle_genes,
    pig %in% pig_genes,
    sheep %in% sheep_genes,
    dog %in% dog_genes
  )
write.table(filtered_data,file = "5species_neg_BW_top1_50kb_homo_genes.txt",sep = '\t',row.names = F)

homogene <- read.table("~/project/10_RER/result/RER/RERconverge_res_20250708.txt",sep = '\t',header = T)
homogene <- homogene[,c(9,10,11,12,13)]
homogene[homogene == "-"] <- NA
homogene <- na.omit(homogene)
filtered_data <- homogene %>%
  filter(
    chicken %in% chicken_genes,
    cattle %in% cattle_genes,
    pig %in% pig_genes,
    sheep %in% sheep_genes
  )
write.table(filtered_data,file = "4species_neg_BW_top1_50kb_homo_genes.txt",sep = '\t',row.names = F)
