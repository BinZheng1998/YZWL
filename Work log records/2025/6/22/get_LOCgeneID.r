setwd("~/project/01_evolution/convergent_evolution/07_result/20250622_LOCgene_selection/BW/")
chicken <- read.table('../../20250605_DCMS_res/BW/chicken/chicken_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(chicken[[5]], ",")
all_genes <- unlist(genes_list)
unique_genes <- unique(all_genes)
symbol <- gsub(" ", "", as.character(unique_genes))
final_genes <- unique(symbol)
final_genes
filtered_genes <- final_genes[grepl("^STRG\\.|^LOC", final_genes)]
filtered_genes
write.table(filtered_genes,"chicken_BW_pos_top005_50kb_LOCgene.txt",row.names = F,col.names = F)

dog <- read.table('../../20250605_DCMS_res/BW/dog/dog_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(dog[[5]], ",")
all_genes <- unlist(genes_list)
unique_genes <- unique(all_genes)
symbol <- gsub(" ", "", as.character(unique_genes))
final_genes <- unique(symbol)
final_genes
filtered_genes <- final_genes[grepl("^STRG\\.|^LOC", final_genes)]
filtered_genes
write.table(filtered_genes,"dog_BW_pos_top005_50kb_LOCgene.txt",row.names = F,col.names = F)

pig <- read.table('../../20250605_DCMS_res/BW/pig/pig_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(pig[[5]], ",")
all_genes <- unlist(genes_list)
unique_genes <- unique(all_genes)
symbol <- gsub(" ", "", as.character(unique_genes))
final_genes <- unique(symbol)
final_genes
filtered_genes <- final_genes[grepl("^STRG\\.|^LOC", final_genes)]
filtered_genes
write.table(filtered_genes,"pig_BW_pos_top005_50kb_LOCgene.txt",row.names = F,col.names = F)

cattle <- read.table('../../20250605_DCMS_res/BW/cattle/cattle_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(cattle[[5]], ",")
all_genes <- unlist(genes_list)
unique_genes <- unique(all_genes)
symbol <- gsub(" ", "", as.character(unique_genes))
final_genes <- unique(symbol)
final_genes
filtered_genes <- final_genes[grepl("^STRG\\.|^LOC", final_genes)]
filtered_genes
write.table(filtered_genes,"cattle_BW_pos_top005_50kb_LOCgene.txt",row.names = F,col.names = F)

sheep <- read.table('../../20250605_DCMS_res/BW/sheep/sheep_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(sheep[[5]], ",")
all_genes <- unlist(genes_list)
unique_genes <- unique(all_genes)
symbol <- gsub(" ", "", as.character(unique_genes))
final_genes <- unique(symbol)
final_genes
filtered_genes <- final_genes[grepl("^STRG\\.|^LOC", final_genes)]
filtered_genes
write.table(filtered_genes,"sheep_BW_pos_top005_50kb_LOCgene.txt",row.names = F,col.names = F)
