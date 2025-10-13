import warnings
warnings.filterwarnings('ignore')
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import spateo as st
import pyvista as pv
import anndata as ad
import spateo as st
import os
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import harmonypy as hm
import colorcet as cc
from colorcet.plotting import swatch
os.system('/usr/bin/Xvfb :99 -screen 0 1024x768x24 &')
os.environ['DISPLAY'] = ':99'
%load_ext autoreload
%autoreload 2

def add_exp_label(adata, embryo_pc, marker):
    pc_index=embryo_pc.point_data["obs_index"].tolist()
    expression = adata[pc_index, marker].X.toarray().flatten()
    norm_expression = expression / expression.max()
    ## 此次可修改绘图颜色，使用归一化的RGB值
    low_color = np.array([0.0, 0.0, 1.0])  # 背景：浅灰色RGB
    high_color = np.array([1.0, 0.0, 0.0])  # 表达颜色：纯红色
    rgb_values = low_color + (high_color - low_color) * norm_expression[:, np.newaxis]
    alpha_values = 0.05 + (1 - 0.05) * norm_expression  # 保持alpha在[0.05, 1]范围
    rgba_array = np.hstack([rgb_values, alpha_values.reshape(-1, 1)])
    embryo_pc.point_data["tissue_rgba"] = rgba_array
    return embryo_pc

def add_celltype_label(adata, embryo_pc, target_celltype):
    is_target = adata.obs['mapped_celltype'].values == target_celltype
    rgba = np.empty((len(is_target), 4), dtype=np.float32)
    ## 此次可修改绘图颜色，使用归一化的RGB值，第四列为透明度
    rgba[is_target] = [0.08, 0.17, 0.55, 1.0]
    rgba[~is_target] = [0.7, 0.7, 0.7, 0.3]
    embryo_pc.point_data["tissue_rgba"] = rgba
    return embryo_pc

def plot_3d_pc(adata, pc, mesh, mode, target, output):
    
    ## 读取anndata格式的空间组数据（3D）
    adata = adata
    embryo_pc = pc
    
    ## 读取blender优化后的mesh文件,并修正格式
    mesh_model = mesh
    ## 设置cell ID
    original_ids = np.arange(mesh_model.n_cells, dtype=np.int64)
    mesh_model.cell_data["vtkOriginalCellIds"] = original_ids
    ## 设置Point ID
    original_point_ids = np.arange(mesh_model.n_points, dtype=np.int64)
    mesh_model.point_data["vtkOriginalPointIds"] = original_point_ids
    ## 设置颜色
    rgba_value = [0.8627451, 0.8627451, 0.8627451, 0.2]
    tissue_rgba = np.zeros((mesh_model.n_cells, 4), dtype=np.float32) 
    tissue_rgba[:, :] = rgba_value
    mesh_model.cell_data["tissue_rgba"] = tissue_rgba
    ## 设置Tissue names
    tissue_array = np.full(mesh_model.n_cells, "surface", dtype=object)
    mesh_model.cell_data["tissue"] = tissue_array
    if mode == 'celltype':
        target_pc = add_celltype_label(adata, embryo_pc, target)
    elif mode == 'marker':
        target_pc = add_exp_label(adata, embryo_pc, target)
    else:
        print('mode error')
    st.pl.three_d_plot(
        model=st.tdr.collect_models([mesh_model, target_pc]), 
        key="tissue", 
        model_style=["surface", "points"], 
        opacity=[0.2, 1],
        model_size=4,
        colormap=None,
        show_legend=False,
        show_axes=True,
        jupyter="static",
        cpo="xy",
        window_size=(1024, 1024),
        background='white',
        plotter_filename = output
    )


## 读取anndata文件
adata = sc.read('/home/data2025/new.chicken_embryo_E3.5_full_demo.h5ad')
lung_cell = pd.read_table('/home/bzheng/project/03_3D_spatial/04_cell_annotation/new_data_2509/lung/E35_lung.txt')
lung_cell = lung_cell['x'].tolist()
sub_ad = adata[lung_cell, :].copy()
sub_ad
sub_ad.obs['type'] = 'lung'
#sub_ad.obs.loc[sub_ad.obs.index.isin(m8), 'type'] = 'm8'
sub_ad.obs['type'].unique()
sub_ad.obs['type'].value_counts()
color_list = ["#E5E5E5","#0000FF"] ##E5E5E5 为grey90
embryo_pc, plot_cmap = st.tdr.construct_pc(adata=sub_ad.copy(), spatial_key="spatial", groupby="type", 
                                           key_added="tissue", 
                                           colormap=color_list)

## 读取点云文件
#embryo_pc = pv.read('/home/data2025/E3.5_embryo_mesh_model.vtk')
## 读取mesh文件
mesh_model = pv.read('/home/data2025/E3.5_embryo_mesh_model.vtk')

## 当绘制细胞类型分布时（细胞类型应是adata.obs['mapped_celltype']中的其中一个）
#plot_3d_pc(adata, embryo_pc, mesh_model, 'celltype', 'fibroblast', 'E4.5_fibroblast.png')

## 当绘制基因表达分布时（基因应是adata.var_names中的其中一个）
plot_3d_pc(sub_ad, embryo_pc, mesh_model, 'marker', 'FGFR2','/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_FGFR2.html')


#

df = sub_ad.obs.z
df.to_csv('/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_z.csv', index=False, encoding='utf-8')
filtered_ad = sub_ad[sub_ad.obs['z'] >= 120].copy()
color_list = ["#E5E5E5","#0000FF"] ##E5E5E5 为grey90
embryo_pc, plot_cmap = st.tdr.construct_pc(adata=filtered_ad.copy(), spatial_key="spatial", groupby="type", 
                                           key_added="tissue", 
                                           colormap=color_list)
mesh_model = pv.read('/home/data2025/E3.5_embryo_mesh_model.vtk')
plot_3d_pc(filtered_ad, embryo_pc, mesh_model, 'marker', 'FGFR2','/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_leftLung.html')
#保存 cell ID
cell_ids = filtered_ad.obs_names
df = pd.DataFrame({'cell_barcode': cell_ids})
df.to_csv('/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_leftLung_cellID.csv',encoding='utf-8')

filtered_ad = sub_ad[sub_ad.obs['z'] <= 106].copy()
color_list = ["#E5E5E5","#0000FF"] ##E5E5E5 为grey90
embryo_pc, plot_cmap = st.tdr.construct_pc(adata=filtered_ad.copy(), spatial_key="spatial", groupby="type", 
                                           key_added="tissue", 
                                           colormap=color_list)
mesh_model = pv.read('/home/data2025/E3.5_embryo_mesh_model.vtk')
plot_3d_pc(filtered_ad, embryo_pc, mesh_model, 'marker', 'FGFR2','/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_rightLung.html')
#保存 cell ID
cell_ids = filtered_ad.obs_names
df = pd.DataFrame({'cell_barcode': cell_ids})
df.to_csv('/home/bzheng/project/03_3D_spatial/02_result/202510-res/E35_rightLung_cellID.csv',encoding='utf-8')
