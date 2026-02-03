setwd('~/project/03_3D_spatial/02_result/260107_E25_to_E25_mesenchyme_monocle3/')
library(Seurat)
library(ggunchull)
library(dplyr)
options(future.globals.maxSize = 50 * 1024^3)
df1 <- readRDS('../../00_data/lung/E25_lung_mesenchymal_1216.rds')
df1$celltype <- 'E25 Lung Mesenchyme'
df1$Stage <- 'E25'
head(df1)
FeaturePlot(df1,features = 'BMP4')
DimPlot(df1,group.by = 'lung_celltype',label = T,raster=FALSE)+NoLegend()
df1 <- subset(df1,subset = lung_celltype == "WNT2+ MSCs")

df2 <- readRDS('../../00_data/lung/E35_lung_mesechymal_251220.rds')
df2$celltype <- 'E35 Lung Mesenchyme'
df2$Stage <- 'E35'
df2 <- NormalizeData(df2,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.1,0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(df2,group.by = 'Spatial_snn_res.0.5',label = T,raster=FALSE)+NoLegend()
FeaturePlot(df1,features = 'NOTCH1')
df2$celltype[df2$Spatial_snn_res.0.5 %in% c('3','4')] <- 'Proximal Mesenchyme'
df2$celltype[df2$Spatial_snn_res.0.5 %in% c('0','1','2')] <- 'Distal Mesenchyme'
DimPlot(df2,group.by = 'celltype',label = T,raster=FALSE)+NoLegend()

head(df2)

obj <- merge(x = df1, 
            y = df2, 
            add.cell.ids = c("E25_mes","E35_mes"), 
            project = "Chicken_Lung_Dev")
DefaultAssay(obj) <- 'Spatial'
obj <- NormalizeData(obj,assay = 'Spatial') %>%
  FindVariableFeatures(nfeatures = 2000,selection.method = 'vst') %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims=1:30) %>%  # 减少维度
  FindClusters(resolution=c(0.5,1)) %>%  # 使用Leiden
  RunUMAP(dims=1:30)
DimPlot(obj,group.by = 'celltype')+ggtitle('k.anchor 1 k.weight 10')

obj_sub<-obj
#obj_sub<- UpdateSeuratObject(obj_sub)
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = RPCAIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.rpca",k.anchor= 2,k.weight=10,
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.rpca", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.rpca', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub,group.by = 'celltype')+ggtitle('k.anchor 1 k.weight 10')

#harmony
obj_sub[["Spatial"]] <- split(obj_sub[["Spatial"]], f = obj_sub$Stage)
obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub <- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:30, reduction = "pca")
obj_sub <- IntegrateLayers(object = obj_sub, method = HarmonyIntegration, orig.reduction = "pca", 
                       new.reduction = "integrated.harmony",theta = 0,lambda = 5,
                       verbose = FALSE)
obj_sub[["Spatial"]] <- JoinLayers(obj_sub[["Spatial"]])
obj_sub <- FindNeighbors(obj_sub, reduction = "integrated.harmony", dims = 1:30)
obj_sub <- FindClusters(obj_sub, resolution = 1)
obj_sub <- RunUMAP(obj_sub, reduction = 'integrated.harmony', dims = 1:30, return.model = T, verbose = F)
DimPlot(obj_sub,group.by = 'celltype')+ggtitle('theta = 0.75,lambda = 1')

library(ggsci)
library(scales)
show_col(pal_npg("nrc")(10))
my_colors <- c(
  "E25 Lung Mesenchyme"  = "#E64b35ff",  # 红色
  "Distal Mesenchyme"  = "#4dbbd5ff",  # 蓝色
  "Proximal Mesenchyme" = "#00a087ff"   # 绿色
)
p1<-DimPlot(obj_sub,group.by = 'celltype',label = T,label.size = 5,cols = my_colors)+ 
    labs(x = '',y='UMAP2',title = '')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(size=15))+
  NoLegend()
p1
FeaturePlot(obj_sub,features = "ADCY8")


###monocle3
#monocle3
library(monocle3)
data <- GetAssayData(obj_sub,assay = 'Spatial',slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('cds.umap')


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(obj_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="Spatial_snn_res.1") + ggtitle('int.umap')


cds <- cluster_cells(cds,cluster_method = 'louvain')
cds <- learn_graph(cds)
p<- plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)
p
p + geom_vline(xintercept = seq(-2.5,-1,0.5)) + geom_hline(yintercept = seq(-4.5,-3.5,0.5))
embed <- data.frame(Embeddings(obj_sub, reduction = "umap"))
embed <- subset(embed, umap_1 > -2.5 & umap_1 < -1 & umap_2 > -4.5 & umap_2 < -3.5)
root.cell <- rownames(embed)
cds <- order_cells(cds, root_cells = root.cell)
#pdf('E25_E35_lung_epithelium_monocle3_pseudotime_umap.pdf',width = 7, height = 5)
p2<-plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 2,
          trajectory_graph_segment_size = 1)+ 
    #labs(x = 'UMAP1',y='UMAP2')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank())+
  NoLegend()
p2
#dev.off()

p1+p2

cds_sub1 <- choose_graph_segments(cds,clear_cds = F,return_list = T)
cds_sub1 <- cds[,cds_sub1$cells]

cds_sub2 <- choose_graph_segments(cds,clear_cds = F,return_list = T)
cds_sub2 <- cds[,cds_sub2$cells]


p3<-plot_cells(cds_sub1, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 1,
          trajectory_graph_segment_size = 1)+ 
  labs(x = 'UMAP1',y='UMAP2')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.title = element_text(size=15))+
  NoLegend()
p3

p4<-plot_cells(cds_sub2, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
           label_leaves = FALSE,  label_branch_points = FALSE,cell_size = 1,
           group_label_size = 1,
          trajectory_graph_segment_size = 1)+ 
  labs(x = 'UMAP1',y='')+
  theme_classic()+
  theme(axis.ticks = element_blank(),
        legend.title = element_text(size = 15),
        #legend.position = 
        legend.text = element_text(size = 15),
        axis.text = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_text(size=15))+
  NoLegend()
p4

#pdf('E25_E35_lung_epithelium_monocle3_pseudotime_umap.pdf',width = 10, height = 8)
p1+p2+p3+p4
#dev.off()



Track_genes1 <- graph_test(cds_sub1, neighbor_graph="principal_graph", cores=20)
#write.table(Track_genes,'E45_lung_mesenchymal_monocle3_res.txt',sep = '\t')
Track_genes_sig1 <- Track_genes1 %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig1 <- c('ADCY8','MLLT3','HOXB6','LRP1B')
p5<-plot_genes_in_pseudotime(cds_sub1[Track_genes_sig1,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
    labs(x='')+
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        strip.text = element_text(face = "italic", size = 12),
        axis.title = element_text(size=15))+NoLegend()
p5

Track_genes2 <- graph_test(cds_sub2, neighbor_graph="principal_graph", cores=20)
#write.table(Track_genes,'E45_lung_mesenchymal_monocle3_res.txt',sep = '\t')
Track_genes_sig2 <- Track_genes2 %>% top_n(n=30, morans_I) %>%
  pull(gene_short_name) %>% as.character()
Track_genes_sig2 <- c('HAND2','AKR1D1','NUAK1','WNT5A')
p6<-plot_genes_in_pseudotime(cds_sub2[Track_genes_sig2,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
  theme_classic()+
  labs(x='Pseudotime')+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        strip.text = element_text(face = "italic", size = 12),
        axis.title = element_text(size=15))+NoLegend()
p6

p <- p1+p2+p3+p4
p

library(patchwork)
pdf('E25_E35_lung_mesenchyme_monocle3_pseudotime_umap.pdf',width = 10, height = 14)
p / (p5 / p6) + plot_layout(heights = c(1.75, 1.25))
dev.off()



Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=20)
#write.table(Track_genes,'E45_lung_mesenchymal_monocle3_res.txt',sep = '\t')
Track_genes_sig1 <- c('HBA1','HBBR','HBZ','HBE','HBAD')
pdf('E25_E35_HBZ_HBAD_HBE_HBBR_HBA1_monocle3_pseudotime_umap.pdf',width = 10, height = 3)
p7<-plot_genes_in_pseudotime(cds_sub1[Track_genes_sig1,],color_cells_by = "pseudotime",min_expr=0.5, ncol = 5)+ 
    labs(x='')+
  theme_classic()+
  theme(
        #axis.line = element_line(linewidth = 1),
        axis.text = element_text(size=12),
        strip.text = element_text(face = "italic", size = 12),
        axis.title = element_text(size=15))+NoLegend()
p7
dev.off()

####GO
library(clusterProfiler)
library(org.Gg.eg.db)

genes1 <- Track_genes1 %>% top_n(n=50, morans_I) %>%
  pull(gene_short_name) %>% as.character()
eg <- bitr(genes1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
id <- as.character(eg[, 2])
ego <- enrichGO(gene = id, OrgDb = "org.Gg.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
ego1 <- as.data.frame(ego)
ego1 <- ego1[ego1$pvalue < 0.05,]
top20_pvalue1 <- ego1 %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(20)
head(top20_pvalue1)
#top20_pvalue$Description <- str_to_title(top20_pvalue$Description)

pdf('E25_to_Distal_monocle3_GO_top20.pdf',width = 7,height = 8)
p_1<-ggplot(top20_pvalue1, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)), 
           fill = RichFactor)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#E8ECF1",  # 这是一个很浅的蓝灰色，比纯白更有质感
                      high = "#4dbbd5ff") + 
  
  geom_text(aes(label = Description, x = 0), 
            hjust = 0, nudge_x = 0.1, size = 5, 
            color = "black") + #以此深蓝为底色时，文字建议改白色，或者根据实际情况改黑色
            
  theme_classic() +
  labs(y = '', x = expression(-log[10](italic(P)-value))) +
  theme(axis.text = element_text(size=12),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
p_1
dev.off()


genes2 <- Track_genes2 %>% top_n(n=50, morans_I) %>%
  pull(gene_short_name) %>% as.character()
eg <- bitr(genes2, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Gg.eg.db")
id <- as.character(eg[, 2])
ego <- enrichGO(gene = id, OrgDb = "org.Gg.eg.db", ont = "all",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                readable = TRUE)
ego2 <- as.data.frame(ego)
ego2 <- ego2[ego2$pvalue < 0.05,]

top20_pvalue2 <- ego2 %>%
  arrange(pvalue) %>%  # 按pvalue升序排列（最小的pvalue最显著）
  head(20)
head(top20_pvalue2)
#top20_pvalue$Description <- str_to_title(top20_pvalue$Description)

pdf('E25_to_Proximal_monocle3_GO_top20.pdf',width = 7,height = 8)
p_2<- ggplot(top20_pvalue2, 
       aes(x = -log10(pvalue),
           y = reorder(Description, -log10(pvalue)), 
           fill = RichFactor)) +
  geom_bar(stat = 'identity', width = 0.8) +
  scale_x_continuous(expand = c(0,0))+
  scale_fill_gradient(low = "#E8ECF1",  # 这是一个很浅的蓝灰色，比纯白更有质感
                      high = "#00a087ff") + 
  
  geom_text(aes(label = Description, x = 0), 
            hjust = 0, nudge_x = 0.1, size = 5, 
            color = "black") + #以此深蓝为底色时，文字建议改白色，或者根据实际情况改黑色
            
  theme_classic() +
  labs(y = '', x = expression(-log[10](italic(P)-value))) +
  theme(axis.text = element_text(size=12),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
p_2
dev.off()

pdf('E25_to_E35_monocle3_GO_top20.pdf',width = 14,height = 8)
p_1+p_2
dev.off()
