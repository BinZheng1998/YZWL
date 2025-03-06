library(Seurat)
library(future)
library(future.apply)
library(networkD3)
library(ggplot2)
library(dplyr)
library(ggalluvial)

subsample <- function(obj){
    cluster_lst <- unique(obj$seurat_clusters)
    sub_cell <- c()
    for (i in 1:length(cluster_lst)){
        obj_sub <- subset(obj, seurat_clusters == cluster_lst[i])
        if(length(Cells(obj_sub)) > 300){
            obj_sub <- obj_sub[, sample(1:ncol(obj_sub),round(ncol(obj_sub)/10))]
            select_cell <- Cells(obj_sub)
        }else{
            select_cell <- Cells(obj_sub)
        }
        sub_cell <- c(sub_cell, select_cell)
    }
    obj <- obj[,sub_cell]
    return(obj)
}


do_Integrate <- function(obj_1, obj_2, prefix){
    seurat_list <- list(obj_1 = obj_1, obj_2 = obj_2)
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA'
    x <- NormalizeData(x, verbose = FALSE)  #数据归一化
    x <- FindVariableFeatures(object = x, verbose = FALSE)
    })
    features <- SelectIntegrationFeatures(object.list = seurat_list)

    #对两个Seurat对象进行标准化和PCA降维
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    return(x)
    })

    # 找到两个数据集的整合锚点
    anchors <- FindIntegrationAnchors(object.list = seurat_list, reduction = "rpca",dims = 1:50)
    obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

    #对整合后的数据进行标准化、PCA、UMAP降维
    obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
    obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
    obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, n.components = 3, min.dist = 0.75)

    emb = data.frame(Embeddings(object = obj.integrated, reduction = "umap"))
    saveRDS(emb, file=paste0("/home/bzheng/project/03_3D_spatial/02_result/heart/20250306-lineage/result/",prefix,".Integrate_umap3.rds"))

    #为两个Seurat对象添加群体和时间标签
    obj_1$subcluster <- paste0("E3.5-", obj_1$subcluster)
    obj_1$Anno = as.vector(obj_1$subcluster)
    obj_1$day = "pre"
    anno1 = obj_1[[]][,c("day", "Anno")]

    obj_2$subcluster <- paste0("E4.5-", obj_2$subcluster)
    obj_2$Anno = as.vector(obj_2$subcluster)
    obj_2$day = "nex"
    anno2 = obj_2[[]][,c("day", "Anno")]

    anno = rbind(anno1, anno2)  # 合并两个时间点的注释信息
    if(nrow(emb) != nrow(anno)){
        print("Error!")
        print(xxx)
    }
    anno = anno[rownames(emb),]  # 按行名排序注释
    res = createLineage_Knn(emb, anno,  k_neigh = 5)
    saveRDS(res, file=paste0("/home/bzheng/project/03_3D_spatial/02_result/heart/20250306-lineage/result/",prefix,".Integrate_Knn_umap_permutation.rds"))
}


createLineage_Knn <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5){
    library(FNN)
    print(dim(emb))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
    pd$state = pd$Anno  # 将群体信息保存在`state`中

    res = list()

    rep_i = 1

    while(rep_i < (replication_times+1)){

        sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))

        emb_sub = emb[sampling_index,] # 获取子集的UMAP嵌入
        pd_sub = pd[sampling_index,] # 获取子集的注释

        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",] # 选取“pre”时间点的细胞
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",] # 选取“nex”时间点的细胞
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]

        pre_state_min = min(table(as.vector(pd_sub1$state))) # 获取最小群体大小

        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min                                                                     
        print(k_neigh)
        }

        if (pre_state_min < 3){
            next
        }

        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index

        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
            tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))

        # 计算邻居之间的状态转化概率
        tmp2 <- matrix(NA,length(state2),length(state1))
        for(i in 1:length(state2)){
            x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
            for(j in 1:length(state1)){                                                                     
                    tmp2[i,j] <- sum(x==state1[j])
            }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1

        res[[rep_i]] = tmp2  # 将结果存储到res列表中

        rep_i = rep_i + 1

    }
    return(res)
}


time_point = paste0("E", c(2.5, 3.5, 4.5))
kk = 1
time_1 = time_point[kk]
time_2 = time_point[kk + 1]
time_3 = time_point[kk + 2]

obj_1 <- readRDS("/home/bzheng/project/03_3D_spatial/02_result/heart/data/SC/E2.5.SC_C27_heart_new.rds")
obj_2 <- readRDS("/home/bzheng/project/03_3D_spatial/02_result/heart/data/SC/E3.5.SC_C22_C31_heart_new.rds")
obj_3 <- readRDS("/home/bzheng/project/03_3D_spatial/02_result/heart/data/SC/E4.5.SC_C33_heart_new.rds")

#obj_1 <- subsample(obj_1)
#obj_2 <- subsample(obj_2)
#obj_3 <- subsample(obj_3)

obj_1@meta.data$group <- paste0(time_1, "_", obj_1@meta.data$orig.ident)
obj_2@meta.data$group <- paste0(time_2, "_", obj_2@meta.data$orig.ident)
obj_3@meta.data$group <- paste0(time_3, "_", obj_3@meta.data$orig.ident)

do_Integrate(obj_1, obj_2, prefix = 'E2.5vsE3.5')
do_Integrate(obj_2, obj_3, prefix = 'E3.5vsE4.5')
