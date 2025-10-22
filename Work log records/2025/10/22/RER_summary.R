setwd('~/project/01_evolution/convergent_evolution/07_result/20251016_5spcies_analysis/7.result//')
library(ggplot2)
df <- read.table('RER_perms_res_251019.txt',sep = '\t',header = T)
df$human_gene <- rownames(df)
df$human_gene <- gsub('_trimmed','',df$human_gene)

orthology_gene <- read.table('../2.RBBH/rbh_combined.txt',sep = '\t',header = T)
df1 <- merge(df,orthology_gene,by.x="human_gene",by.y="human")
df1$human <- df1$human_gene
gene <- read.table('../1.protein_data/10species_info.txt',sep = '\t',header = T)

cols_to_replace <- c(1, 9:18)
for (col in cols_to_replace) {
  match_index <- match(df1[[col]], gene$protein_id1)
  df1[[col]] <- ifelse(is.na(match_index), 
                       df1[[col]], 
                       gene$gene_id[match_index])
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Gg.eg.db)


######################RER acceleration
acc_res <- df1[df1$P < 0.05 & df1$Rho > 0 ,]
genes <- unique(acc_res$human_gene)
genes
eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
id <- as.character(eg[, 2])
ego <- enrichGO(gene = id, OrgDb = "org.Hs.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
#ego_simplified <- simplify(ego, cutoff =0.7, by ="pvalue", select_fun = min)
ego <- as.data.frame(ego)
ego$Description <- str_to_title(ego$Description)
ego <- ego[ego$pvalue < 0.05, ]
top25_pvalue <- ego %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(25)
head(top25_pvalue)
top25_pvalue$Description <- str_to_title(top25_pvalue$Description)

pdf('RER_acceleration_GO_top25.pdf',width = 10,height = 7)
ggplot(top25_pvalue, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)),  # 按p-value排序
           fill = RichFactor)) +
  scale_fill_distiller(palette = 'RdPu', direction = 1) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_classic() +
  labs(y = '', x = '-log10(p-value)')+
  theme(axis.text = element_text(size=12))
dev.off()

######################RER deceleration
del_res <- df1[df1$P < 0.05 & df1$Rho < 0 ,]
genes <- unique(del_res$human_gene)
genes
eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
id <- as.character(eg[, 2])

ego <- enrichGO(gene = id, OrgDb = "org.Hs.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
#ego_simplified <- simplify(ego, cutoff =0.7, by ="pvalue", select_fun = min)
ego <- as.data.frame(ego)
ego$Description <- str_to_title(ego$Description)
ego <- ego[ego$pvalue < 0.05, ]
top25_pvalue <- ego %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(25)
head(top25_pvalue)
top25_pvalue$Description <- str_to_title(top25_pvalue$Description)
pdf('RER_deceleration_GO_top25.pdf',width = 10,height = 7)
ggplot(top25_pvalue, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)),  # 按p-value排序
           fill = RichFactor)) +
  scale_fill_distiller(palette = 'RdPu', direction = 1) +
  geom_bar(stat = 'identity', width = 0.8) +
  theme_classic() +
  labs(y = '', x = '-log10(p-value)')+
  theme(axis.text = element_text(size=12))
dev.off()
