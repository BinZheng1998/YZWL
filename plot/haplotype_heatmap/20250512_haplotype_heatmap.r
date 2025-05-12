setwd("E:/02-群体进化/07-结果/20250429-单倍型热图/pig/")
Group <- read.table("pig_BW_metadata.txt",stringsAsFactors=F)
names(Group) <- c("Sample","Group")
Data <- read.table("pig_BW_FANCA.txt",header = T,stringsAsFactors=F)
df <- merge(Group,Data,by="Sample")
df <- df[order(df[,2]),]

group <- df[,c(1,2)]
group <- group[,-1]
group <- data.frame(group)
rownames(group)<-df$Sample
names(group) <- "Group"

data=df[,-2]
rownames(data) <- data$Sample
data=data[,-1]
data=data.matrix(data)
an_row <- data.frame(Group=as.vector(group$Group))
rownames(an_row)=rownames(data)
ra_col<-list(Group=c(BW_high="#18A2CA",BW_low="#E20593"))
library(pheatmap)
library(ComplexHeatmap)
pheatmap(data,cluster_row=F,cluster_col=F,legend=F,
         show_rownames=FALSE,show_colnames=FALSE,
         annotation_legend=T,annotation_names_row=F,
         annotation_row=an_row,annotation_colors=ra_col,
         gaps_row=486,
         color=c("grey","#BF3826","#F7F8D5"))

genotype <- sort(unique(as.character(data)))
colors <- c("grey","#BF3826","#F7F8D5")
names(colors) <- genotype
head(colors)

pdf("pig_FANCA.pdf",width = 8,height = 6)
res <-Heatmap(
  data,
  border = FALSE,
  col = colors,
  name = " ",
  heatmap_legend_param = list(color_bar="discrete",
                              at=c(0,1,2),nrow=1,
                              #direction="horizontal",
                              labels=c("Homozygous reference","Heterozygous variant","Homozygous variant")),
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = rowAnnotation(
    group = anno_block(gp = gpar(fill = c("grey95", "grey95"),col=NA), 
                       labels = c("BW high", "BW low")             
    )
  ),
  row_split = df$Group,  # 确保 df$Group 是因子且水平与颜色顺序匹配
  row_title = NULL,
  show_heatmap_legend = T
)
draw(res,heatmap_legend_side="bottom")
dev.off()
