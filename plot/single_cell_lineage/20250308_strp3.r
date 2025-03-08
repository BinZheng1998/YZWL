library(networkD3)
library(ggplot2)
library(dplyr)
library(ggalluvial)

edge_prob <- read.table("/home/bzheng/project/03_3D_spatial/02_result/heart/20250306-lineage/result/SP/E2.5_VS_E3.5_VS_E4.5.SP.edge_prob.txt", header = FALSE, sep = "\t", col.names = c("from", "to", "weight"))
nodes <- data.frame(
  name=c(as.character(edge_prob$from), 
         as.character(edge_prob$to)) %>% unique()
)
edge_prob$IDsource <- match(edge_prob$from, nodes$name)-1 
edge_prob$IDtarget <- match(edge_prob$to, nodes$name)-1
p <- sankeyNetwork(Links = edge_prob, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "weight", NodeID = "name", 
                   sinksRight=FALSE)

saveNetwork(p,"/home/bzheng/project/03_3D_spatial/02_result/heart/20250306-lineage/result/SP/20250308_SP_Heart_E2.5_vs_E3.5_vs_E4.5.subcluster_lineage.html")
