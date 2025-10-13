import anndata as ad
import spateo as st
import os
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
import harmonypy as hm
import colorcet as cc
from colorcet.plotting import swatch
os.system('/usr/bin/Xvfb :99 -screen 0 1024x768x24 &')
os.environ['DISPLAY'] = ':99'
%load_ext autoreload
%autoreload 2

data = sc.read('/home/data2025/new.chicken_embryo_E3.5_full_demo.h5ad')
lung_cell = pd.read_table('/home/bzheng/project/03_3D_spatial/04_cell_annotation/new_data_2509/lung/E35_lung.txt')
lung_cell = lung_cell['x'].tolist()
m19 = pd.read_table('/home/bzheng/project/03_3D_spatial/04_cell_annotation/new_data_2509/lung/E35_lung_M8.txt')
m19 = m19['x'].tolist()
sub_ad = data[lung_cell, :].copy()
sub_ad
sub_ad.obs['type'] = 'lung'
sub_ad.obs.loc[sub_ad.obs.index.isin(m19), 'type'] = 'm8'
sub_ad.obs['type'].unique()
sub_ad.obs['type'].value_counts()
color_list = ["#808080","#0000FF"]
embryo_pc, plot_cmap = st.tdr.construct_pc(adata=sub_ad.copy(), spatial_key="spatial", groupby="type", key_added="tissue", colormap=color_list)

import pyvista as pv
mesh_model = pv.read('/home/data2025/E3.5_embryo_mesh_model.vtk')
## 设置cell ID
original_ids = np.arange(mesh_model.n_cells, dtype=np.int64)
mesh_model.cell_data["vtkOriginalCellIds"] = original_ids
## 设置Point ID
original_point_ids = np.arange(mesh_model.n_points, dtype=np.int64)
mesh_model.point_data["vtkOriginalPointIds"] = original_point_ids
## 设置颜色
rgba_value = [0.8627451, 0.8627451, 0.8627451, 0.2]
tissue_rgba = np.zeros((mesh_model.n_cells, 4), dtype=np.float32)  # 使用 float32 节省内存
tissue_rgba[:, :] = rgba_value
mesh_model.cell_data["tissue_rgba"] = tissue_rgba
## 设置Tissue names
tissue_array = np.full(mesh_model.n_cells, "surface", dtype=object)
mesh_model.cell_data["tissue"] = tissue_array

st.pl.three_d_plot(
    model=st.tdr.collect_models([mesh_model, embryo_pc]), 
    key="tissue", 
    model_style=["surface", "points"],
    opacity=[0.2, 1],
    model_size=2,
    colormap=None,
    show_legend=True,
    show_axes=True,
    jupyter="static",
    cpo="xy",
    window_size=(1024, 1024),
    background='white',
    plotter_filename = f'/home/bzheng/project/03_3D_spatial/02_result/202510-res/E3.5.lung_type.embryo_pc_model.html'
    #filename = f'E4.5.celltype.embryo_pc_model.pdf'
)
