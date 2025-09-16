library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
data <- read.table('celltype_tissue.txt',sep = '\t',header = T)
data$colors <- paste0('#',data$colors)

# Assuming data$celltype_gut has 35 unique labels
n_labels <- length(unique(data$celltype_gut))
stopifnot(n_labels == 35) # Verify the number of labels


# 35种颜色的16进制代码
colors <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", 
  "#8C564B", "#E377C2", "#7F7F5F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5",
  "#393B79", "#637939", "#8C6D31", "#843C39", "#7B4173",
  "#5254A3", "#6BAED6", "#9E9AC8", "#969696", "#6B6ECF",
  "#FD8D3C", "#E6550D", "#31A399", "#756BB1", "#636383"
)
colors1<-data$colors
numbered_labels <- data$celltype_gut

library(ComplexHeatmap)
library(grid)


lgd <- Legend(
  labels = numbered_labels, 
  type = "points", 
  pch = 21,
  legend_gp = gpar(col = colors1, fill = colors1),
  ncol = 4, 
  row_gap = unit(2, "mm"),
  title_position = "topcenter",
  background = "white",
  size = unit(5, "mm"), # 圆的大小
  labels_gp = gpar(fontsize = 10), # 标签字体大小
  grid_height = unit(6, "mm"), # 图例项高度
  grid_width = unit(4, "mm"), # 减小宽度以拉近标签
  graphics = lapply(seq_along(numbered_labels), function(i) {
    function(x, y, w, h) {
      grid.circle(x = x - unit(1, "mm"), y = y, r = unit(2.5, "mm"), 
                  gp = gpar(col = colors1[i], fill = colors1[i]))
      grid.text(
        label = as.character(i - 1), # 从 0 开始，与标签对应
        x = x - unit(1, "mm"), # 减小偏移量
        y = y,
        gp = gpar(fontsize = 8, col = "white",fontface = "bold")
      )
    }
  })
)
pdf('Gut_legend.pdf',width = 12,height = 4)
draw(lgd)
dev.off()
