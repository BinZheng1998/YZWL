setwd("E:/02-群体进化/07-结果/20250603-富集分析聚类/")
library(simplifyEnrichment)
library(org.Gg.eg.db)
library(org.Hs.eg.db)
BiocManager::install("org.Hs.eg.db")
library(ComplexHeatmap)
ego <- read.table("E:/02-群体进化/07-结果/202504-富集分析/GO-KEGG/BW/chicken/chicken_BW_pos_DCMS_top5_merged_regions_50kb_GO.txt",header = T)

go_id_bp <- ego[ego$ONTOLOGY == "BP", "ID"]
go_id_cc <- ego[ego$ONTOLOGY == "CC", "ID"]
go_id_mf <- ego[ego$ONTOLOGY == "MF", "ID"]
length(go_id_bp);length(go_id_cc);length(go_id_mf)
head(go_id_bp)
mat <- GO_similarity(go_id_mf,ont = "MF",db="org.Gg.eg.db")
pdf("heatmap_chicken_pos_BW_top_50kb_MF.pdf", width = 12, height = 8)
simplifyGO(mat, plot = T)
dev.off()
