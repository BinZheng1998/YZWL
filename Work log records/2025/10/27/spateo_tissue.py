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
import pyvista as pv
from colorcet.plotting import swatch
os.system('/usr/bin/Xvfb :99 -screen 0 1024x768x24 &')
os.environ['DISPLAY'] = ':99'
%load_ext autoreload
%autoreload 2

adata = sc.read('/home/data2025/new.chicken_embryo_E4.5_full_demo.h5ad')
mesh_model = pv.read('/home/data2025/E4.5_embryo_mesh_model.vtk')

# 设置网格的细胞和点 ID
mesh_model.cell_data["vtkOriginalCellIds"] = np.arange(mesh_model.n_cells, dtype=np.int64)
mesh_model.point_data["vtkOriginalPointIds"] = np.arange(mesh_model.n_points, dtype=np.int64)

# 设置网格默认颜色
rgba_value = [0.8627451, 0.8627451, 0.8627451, 0.2]
mesh_model.cell_data["tissue_rgba"] = np.full((mesh_model.n_cells, 4), rgba_value, dtype=np.float32)
mesh_model.cell_data["tissue"] = np.full(mesh_model.n_cells, "surface", dtype=object)

# 颜色映射
color_map = {
    'Unknown': [0.87451, 0.345098, 1.0],
    'Tail Bud': [0.54902, 0.235294, 1.0],
    'Blood Vessels': [0.007843, 0.533333, 0.0],
    'Heart': [0.0, 0.67451, 0.780392],
    'Macrophages': [0.596078, 1.0, 0.0],
    'Craniofacial Region': [1.0, 0.647059, 0.188235],
    'Ear': [0.423529, 0.0, 0.309804],
    'Neural Crest': [1.0, 0.498039, 0.819608],
    'Pharyngeal Arch': [0.0, 0.0, 0.615686],
    'Thyroid': [0.52549, 0.439216, 0.407843],
    'Duodenum': [0.0, 0.286275, 0.258824],
    'Foregut': [0.309804, 0.164706, 0.0],
    'Gizzard': [0.0, 0.992157, 0.811765],
    'Gut': [0.145098, 0.4, 0.635294],
    'Hindgut': [0.584314, 0.705882, 0.478431],
    'Liver': [0.752941, 0.015686, 0.72549],
    'Lung': [0.737255, 0.717647, 1.0],
    'Pancreas': [0.156863, 0.0, 0.254902],
    'Proventriculus': [0.862745, 0.701961, 0.686275],
    'Limb Bud': [0.996078, 0.960784, 0.564706],
    'Surface Ectoderm': [0.313725, 0.270588, 0.356863],
    'Allantois': [0.643137, 0.486275, 0.0],
    'Amnion': [1.0, 0.443137, 0.4],
    'Intermediate Mesoderm': [0.247059, 0.505882, 0.431373],
    'Lateral Plate Mesoderm': [0.509804, 0.0, 0.05098],
    'Splanchnic Mesoderm': [0.639216, 0.482353, 0.701961],
    'Diencephalon': [0.203922, 0.305882, 0.0],
    'Eye': [0.607843, 0.894118, 1.0],
    'Hypothalamus': [0.921569, 0.0, 0.466667],
    'Isthmus': [0.176471, 0.0, 0.039216],
    'Mesencephalon': [0.368627, 0.564706, 1.0],
    'Neural Tube': [0.0, 0.780392, 0.12549],
    'Olfactory Bulb': [0.345098, 0.003922, 0.666667],
    'Pineal Gland': [0.0, 0.117647, 0.0],
    'Primitive Streak': [0.603922, 0.278431, 0.0],
    "Rathke's Pouch": [0.588235, 0.623529, 0.65098],
    'Rhombencephalon': [0.607843, 0.258824, 0.360784],
    'Telencephalon': [0.0, 0.121569, 0.196078],
    'Notochord': [1.0, 0.815686, 1.0],
    'Somites': [0.843137, 0.0, 0.0],
    'Cloaca': [0.0, 0.745098, 0.603922],
    'Genital Tubercle': [0.215686, 0.082353, 1.0],
    'Gonad': [0.176471, 0.145098, 0.145098],
    'Mesonephros': [0.784314, 0.768627, 0.0]
}

for tissue, color in color_map.items():
    # 子集 AnnData
    sub_ad = adata[adata.obs['Tissue'] == tissue]
    if sub_ad.n_obs == 0:
        print(f"组织 {tissue} 无细胞，跳过...")
        continue
    
    sub_ad.obs['type'] = tissue
    print(f"组织：{tissue}，细胞数：{sub_ad.n_obs}")
    
    # 构建点云
    embryo_pc, plot_cmap = st.tdr.construct_pc(
        adata=sub_ad.copy(),
        spatial_key="spatial",
        groupby="type",
        key_added="tissue",
        colormap=[color]  # 使用单一颜色列表
    )
    
    # 手动设置点云颜色
    if "tissue" in embryo_pc.point_data:
        tissue_labels = embryo_pc.point_data["tissue"]
        colors = np.array([color_map[tissue] for _ in tissue_labels])
        embryo_pc.point_data["tissue"] = np.hstack([colors, np.ones((colors.shape[0], 1))])  # 使用 "tissue" 存储 RGBA
    else:
        print(f"点云中未找到组织数据：{tissue}")
        continue
    
    # 验证点云
    print(f"组织：{tissue}，点数：{embryo_pc.n_points}")
    print(f"点云组织数据：{embryo_pc.point_data.get('tissue', '无组织数据')[:5]}")
    print(f"点云颜色数据：{embryo_pc.point_data.get('tissue', '无颜色数据')[:5]}")
    
    # 生成 3D 图
    tissue_name = tissue.replace(' ', '_')
    st.pl.three_d_plot(
        model=st.tdr.collect_models([mesh_model, embryo_pc]),
        key="tissue",  # 使用 "tissue" 作为 key
        model_style=["surface", "points"],
        opacity=[0.2, 1],
        model_size=3,
        colormap=None,
        show_legend=False,
        show_axes=True,
        jupyter="static",
        cpo="xy",
        window_size=(1424, 1424),
        background='white',
        plotter_filename=f'/home/bzheng/project/03_3D_spatial/02_result/202510-res/E4.5/eye.subtype.3d.{tissue_name}.html'
    )
    
    print(f"已为 {tissue} 生成图")
