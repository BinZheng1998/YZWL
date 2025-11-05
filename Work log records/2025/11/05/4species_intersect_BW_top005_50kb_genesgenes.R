setwd("~/project/01_evolution/convergent_evolution/07_result/20250801_BW_DCMS/BW_convergent/")
library(dplyr)
animals <- c("chicken", "cattle", "pig", "sheep", "dog")
all_species_genes <- list()
for (animal in animals) {
    message(paste("Processing:", animal, "..."))
    input_file <- paste0("../", animal, "/", animal, "_BW_pos_DCMS_top_005_merged_regions_50kb_genes.txt")
    data <- read.table(input_file, sep = '\t', header = TRUE)
    genes_list <- strsplit(data[[5]], ",")
    all_genes <- unlist(genes_list)
    symbol <- as.character(all_genes)
    symbol <- gsub(" ", "", symbol) # Remove spaces
    final_genes <- sort(unique(symbol), decreasing = FALSE)
    all_species_genes[[animal]] <- final_genes
}


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
#write.table(filtered_data,file = "4species_pos_BW_top5_50kb_homo_genes.txt",sep = '\t',row.names = F)

library(clusterProfiler)
library(org.Hs.eg.db)
eg <- bitr(filtered_data$pig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
id <- as.character(eg[, 2])

ego <- enrichGO(gene = id, OrgDb = "org.Hs.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
ego <- as.data.frame(ego)
ego <- ego[ego$pvalue < 0.05,]

top20_pvalue <- ego %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(20)
head(top20_pvalue)
top20_pvalue$Description <- str_to_title(top20_pvalue$Description)

#top20_pvalue$Description[top20_pvalue$ID == "GO:0045935"] <- "Positive Regulation Of Nucleobase-Containing Compound Metabolic"

pdf('4species_pos_top5_50kb_GO_top20.pdf',width = 10,height = 7)
ggplot(top20_pvalue, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)),  # 按p-value排序
           fill = RichFactor)) +
  scale_fill_distiller(palette = 'RdPu', direction = 1) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_classic() +
  labs(y = '', x = '-log10(p-value)')+
  theme(axis.text = element_text(size=12))
dev.off()


x <- list('Cattle' = cattle_genes,
          'Sheep' = sheep_genes,
          'Chicken' = chicken_genes,
          'Pig' = pig_genes)

ortho_gene <- homogene[,c(1,3,4,5)]
colnames(ortho_gene) <- c('Pig','Sheep','Cattle','Chicken')
ortho_gene[ortho_gene == "None"] <- NA
ortho_gene$OG_ID <- paste0("OG", 1:nrow(ortho_gene))
ortholog_venn_list <- list()
for (species in names(x)) {
  
  # 1. 获取该物种的基因集 (e.g., cattle_genes)
  current_gene_list <- x[[species]]
  
  # 2. 检查直系同源表中是否有该物种的列
  if (!species %in% colnames(ortho_gene)) {
    warning(paste("找不到列:", species, "在直系同源表中, 已跳过"))
    next # 跳过这个物种
  }
  
  # 3. 获取直系同源表中该物种的基因列 (e.g., ortho_gene$Cattle)
  ortholog_column_genes <- ortho_gene[[species]]
  
  # 4. 找到你的基因集与直系同源表的交集, 并获取对应的 OG_ID
  #    (找出 ortho_gene 中, 哪一行的基因 存在于 current_gene_list)
  #    %in% 会自动处理 NA (返回 FALSE)，所以很安全
  matching_ogs <- ortho_gene$OG_ID[ortholog_column_genes %in% current_gene_list]
  
  # 5. 存储唯一的 OG_ID 到新列表中
  #    (使用 unique 以防万一)
  ortholog_venn_list[[species]] <- unique(matching_ogs)
}
print(paste("Cattle OG 数量:", length(ortholog_venn_list$Cattle)))
print(paste("Sheep OG 数量:", length(ortholog_venn_list$Sheep)))
print(paste("Chicken OG 数量:", length(ortholog_venn_list$Chicken)))
print(paste("Pig OG 数量:", length(ortholog_venn_list$Pig)))

library(ggvenn)
ggvenn(ortholog_venn_list,
       show_percentage = FALSE,
       stroke_color = "white",
       fill_color = c("#ffb2b2","#b2e7cb","#b2d4ec","#d3c0e2"),
       set_name_color = c("#ff0000","#4a9b83","#1d6295","#7030a2"))

write.table(filtered_data,'4species_intersect_pos_BW_top5_50kb_homo_genes.txt',row.names = F,sep = '\t')
