library(dbscan)
coords_3d <- E25epi@meta.data[, c("trans_x", "trans_y", "trans_z")]
db_res <- dbscan(coords_3d, eps = 50, minPts = 5)
E25epi$cluster_3d <- db_res$cluster
DimPlot(E25epi,group.by = 'cluster_3d')
source('~/project/03_3D_spatial/03-script/export_clusters_to_ply.R')

export_clusters_to_ply(E25epi,output_dir = './E25epi',coord_cols = c("trans_x", "trans_y", "trans_z"),cluster_col = 'cluster_3d')
