解卷积分析流程

使用cell2location，将单细胞数据的注释映射到空间数据上。具体步骤：
b
1. 数据准备：
    a. 准备anndata格式的单细胞数据和空间组数据，单细胞数据应有注释，空间组数据应有坐标信息（储存在adata.obsm['spatial']）;
    b. 单细胞和空间组数据使用的矩阵（adata.X）应为原始counts矩阵（print(adata.X[1:10,1:10])检查矩阵是否为counts）；

2. 环境配准
    参考工具网站主页说明，配准环境：https://github.com/BayraktarLab/cell2location

3. 准备脚本：
    a. cell2loc运行python源脚本为cell2loc.py，一般不用更改；
    b. 通过shell脚本设置运行cell2loc的参数

        ****/Miniconda3/envs/cell2loc_env/bin/python \   ## cell2loc运行环境
            /share/appspace_data/shared_groups/bgi_liaox_liaoxun/Chicken_embryo/04.cell2loc/cell2loc.py \    ## cell2loc源脚本路径
            -sp /share/home/bgi_liaox/workplace/06.Chicken_embryo/03.regist/plot_exp/anndata/E3.5/D03554C2.bin50.raw.h5ad \     ## 空间组数据
            -ref /share/appspace_data/shared_groups/bgi_liaox_liaoxun/Chicken_embryo/03.SingleCell/00.object/E3.5.SC.C0_C2.leg_bud.recluster.h5ad \    ## 单细胞数据
            -outpfx E3.5-C0_C2_to_D03554C2 \      ## 输出文件前缀
            -label seurat_clusters \              ## 单细胞数据中储存注释或要映射信息的列
            -num 4                                ## 空间组数据中，每一个spot中平均细胞数，bin50设置为4，bin100则设置为16

3. 投递脚本：
    a. 任务应投递到GPU节点中；
