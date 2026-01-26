setwd('~/project/03_3D_spatial/00_data/database/cellchat/')
load('~/project/03_3D_spatial/00_data/database/cellchat/CellChatDB.mouse.rda')
CellChatDB.use <- subsetDB(CellChatDB.mouse)

data <- CellChatDB.use$interaction

orthgene <- read.table('~/project/03_3D_spatial/00_data/database/cellchat/mouse_chicken_ensembl_orthologues.txt',sep = '\t',header = T)
orthgene1 <- orthgene[,c(5,7)]
orthgene1 <- unique(orthgene1)
orthgene1[orthgene1 == ""] <- NA
orthgene2 <- na.omit(orthgene1)
orthgene2 <- unique(orthgene2)
colnames(orthgene2) <- c("Chicken_Gene","Mouse_Gene")
gene_map <- setNames(orthgene2$Chicken_Gene, orthgene2$Mouse_Gene)
trans_gene <- function(genes, map_vec) {
  new_genes <- map_vec[genes]
  return(new_genes)
}

df_inter <- CellChatDB.use$interaction
complex_names <- rownames(CellChatDB.use$complex)
df_inter$ligand <- sapply(df_inter$ligand, function(x) {
  if (x %in% complex_names) return(x) # 如果是复合物名，保持不变
  if (x %in% names(gene_map)) return(gene_map[x]) # 如果有同源基因，替换
  return(NA) # 既不是复合物也没同源基因，标记删除
})

df_inter$receptor <- sapply(df_inter$receptor, function(x) {
  if (x %in% complex_names) return(x)
  if (x %in% names(gene_map)) return(gene_map[x])
  return(NA)
})
df_inter_chicken <- df_inter %>% 
  filter(!is.na(ligand) & !is.na(receptor))
df_inter_chicken$interaction_name_2 <- paste0(df_inter_chicken$ligand, " - ", df_inter_chicken$receptor)
df_inter_chicken$ligand.symbol <- df_inter_chicken$ligand
df_inter_chicken$receptor.symbol <- df_inter_chicken$receptor
df_inter_chicken$ligand.symbol <- gsub("_", ",", df_inter_chicken$ligand.symbol)
df_inter_chicken$receptor.symbol <- gsub("_", ",", df_inter_chicken$receptor.symbol)
df_inter_chicken$interaction_name_2 <- gsub("_", "+", df_inter_chicken$interaction_name_2)

df_complex <- CellChatDB.use$complex
for (col in paste0("subunit_", 1:4)) {
  if (col %in% colnames(df_complex)) {
    df_complex[[col]] <- sapply(df_complex[[col]], function(x) {
      if (x == "" || is.na(x)) return("") # 空值保持空
      if (x %in% names(gene_map)) return(gene_map[x]) # 替换
      return(NA) # 找不到同源基因标记为 NA
    })
  }
}

df_complex_chicken <- na.omit(df_complex)
valid_complexes <- rownames(df_complex_chicken)
valid_entities <- c(orthgene2$Chicken_Gene, valid_complexes)

df_inter_chicken <- df_inter_chicken %>%
  filter(ligand %in% valid_entities & receptor %in% valid_entities)

#修改geneInfo
df_geneInfo <- CellChatDB.use$geneInfo
df_geneInfo$Symbol <- sapply(df_geneInfo$Symbol, function(x) {
  if (x %in% names(gene_map)) return(gene_map[x])
  return(NA)
})
df_geneInfo_chicken <- df_geneInfo[!is.na(df_geneInfo$Symbol), ]
df_geneInfo_chicken <- df_geneInfo_chicken[!duplicated(df_geneInfo_chicken$Symbol), ]

CellChatDB.chicken <- list(
  interaction = df_inter_chicken,
  complex = df_complex_chicken,
  cofactor = CellChatDB.use$cofactor,
  geneInfo = df_geneInfo_chicken
)
save(CellChatDB.chicken, file = "CellChatDB.chicken.rda")
