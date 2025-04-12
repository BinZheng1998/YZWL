import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np

## 01.读取****.denovo.h5ad文件
adata = sc.read('****.denovo.h5ad')
## 02.把解卷积结果写到obs里，方便后面画图
adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
## 03.查看写到obs里的解卷积结果
adata

## 04.画图
outpfx = '/share/appspace_data/shared_groups/bgi_liaox_liaoxun/Chicken_embryo/03.SingleCell/01.result/E3.5-C0_C2_to_Y00560L8'    ## 输出文件路径及前缀
fig, axs = plt.subplots(1, 1, figsize=(15, 15))
sc.pl.spatial(adata, cmap='magma',
                  # 从obs中选择要绘图的列
                  color=[0, 1, 2],
                  ncols=4, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% cd
                  vmin=0, vmax='p99.2',
                  spot_size=42,
                  #spot_size=sizes,
                  show=False
              )
plt.savefig(f"{outpfx}.denovo.cell2loc.pdf")
