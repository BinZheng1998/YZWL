setwd('~/project/01_evolution/convergent_evolution/07_result/20250730_DCMS_GO/BW/chicken/')
library(clusterProfiler)
library(GO.db)
library(org.Gg.eg.db)
pos_gene <- read.table('../../../20250714_DCMS_res/BW/chicken/chicken_BW_pos_DCMS_top_001_merged_regions_50kb_genes.txt',sep = '\t',header = T)
genes_list <- strsplit(pos_gene[[5]], ",")
all_genes <- unlist(genes_list)
symbol <- as.character(all_genes)
symbol <- gsub(" ","",symbol)
genes <- unique(symbol)
genes
eg <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
id <- as.character(eg[, 2])

ego <- enrichGO(gene = id, OrgDb = "org.Gg.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)

ego <- as.data.frame(ego)
ego <- ego[ego$pvalue < 0.05,]
write.table(ego,'chicken_DCMS_pos_top1_50kb_GO.txt',sep = '\t',row.names = F)
library(gground)
library(ggprism)
library(dplyr)
use_pathway <- group_by(ego,ONTOLOGY) %>% top_n(5,wt = -pvalue) %>% group_by(pvalue) %>% head(n=15)%>% top_n(1,wt = Count) %>%
  ungroup() %>% mutate(ONTOLOGY = factor(ONTOLOGY,levels = rev(c("BP","CC","MF")))) %>% arrange(ONTOLOGY,pvalue) %>%
  mutate(Description = factor(Description,levels = Description)) %>% tibble::rowid_to_column('index') %>% 
  mutate(geneID = map_chr(str_split(geneID, "/"), ~ {  #提取前15个基因展示
    selected_genes <- head(.x, 15) 
    str_c(selected_genes, collapse = "/")}))
#width <- 0.25
xaxis_max <- max(-log10(use_pathway$pvalue)) +1
rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  reframe(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -0.7,
    xmax = -0.4,
    ymax = cumsum(n),
    ymin = lag(ymax, default =0) +0.6,
    ymax = ymax +0.4
  )
rect.data
pal <- c('#E64b35b2','#4dbbd5b2','#3c5488b2','#ed6ca4')

p<-ggplot(use_pathway,aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width =0.6, alpha =0.8
  ) +
  geom_text(
    aes(x =0.05, label = Description),
    hjust =0, size =4.5
  ) +
  geom_text(
    aes(x =0.1, label = geneID, colour = ONTOLOGY),
    hjust =0, vjust =3, size =3, fontface ='italic',
    show.legend =FALSE
  ) +
  # 基因数量
  geom_point(
    aes(x = -0.2, size = Count),
    shape =21
  ) +
  geom_text(
    aes(x = -0.2, label = Count)
  ) +
  scale_size_continuous(name ='Count', range = c(5,12)) +
  # 分类标签
  geom_round_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = ONTOLOGY),
    data = rect.data,
    radius = unit(2,'mm'),
    inherit.aes =FALSE
  ) +
  geom_text(
    aes(x = (xmin + xmax) /2, y = (ymin + ymax) /2, label = ONTOLOGY),
    data = rect.data,
    inherit.aes =FALSE
  ) +
  geom_segment(
    aes(x =0, y =0, xend = xaxis_max, yend =0),
    linewidth =1.5,
    inherit.aes =FALSE
  ) +
  labs(y =NULL) +
  scale_fill_manual(name ='Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, xaxis_max,2),
    expand = expansion(c(0,0))
  ) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none'
    #legend.title = element_text()
  )
p
ggsave('chicken_DCMS_pos_top1_50kb_GO.pdf',p,width =10, height =8,dpi=300)

library(aPEAR)
p1<-enrichmentNetwork(ego,colorBy = 'pvalue',
                      colorType = 'pval',
                      nodeSize = 'zScore',
                      fontSize = 4,minClusterSize = 3,
                      drawEllipses = T,
                      verbose = F)
p1
ggsave('chicken_DCMS_pos_top1_50kb_GO_cluster.pdf',p1,width = 16,height = 12,dpi=300)
