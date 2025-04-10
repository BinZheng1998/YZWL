import scanpy as sc
import pandas as pd
import numpy as np
import os
import cell2location
from cell2location.models import RegressionModel
import scvi

def regression_model(adata_ref, 
                    save_path, 
                    batch_key, 
                    labels_key, 
                    accelerator="gpu",
                    sc_model=None,
                    ):
    
    import cell2location
    from cell2location.models import RegressionModel
    import scvi
    if sc_model is None:
        # prepare anndata for the regression model
        cell2location.models.RegressionModel.setup_anndata(
            adata=adata_ref,
            # 10X reaction / sample / batch
            batch_key=batch_key,
            # cell type, covariate used for constructing signatures
            labels_key=labels_key,
            )

        # create and train the regression model
        mod = RegressionModel(adata_ref)
        
        print('the setup data as below')
        mod.view_anndata_setup()

        print('Use all data for training ' 
            '(validation not implemented yet, `train_size=1`)')
        
        mod.train(max_epochs=200, 
                batch_size=2500,
                train_size=1, 
                lr=0.005, 
                accelerator="gpu")

        # plot ELBO loss history during training
        # removing first 20 epochs from the plot
        mod.plot_history(20)

        # Save model
        mod.save(f"{save_path}", overwrite=True)
    
    else:
        mod = cell2location.models.Cell2location.load(f"{save_path}", adata_ref)
        
    # In this section, we export the estimated 
    # cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(adata_ref, 
                                    sample_kwargs={'num_samples': 1000, 
                                                    'batch_size': 2500, 
                                                    'accelerator': 'gpu'}
                                    )
    # Save anndata object with results
    adata_file = f"{save_path}/sc.h5ad"
    adata_ref.write(adata_file)
    print(f'{adata_file} saved!')
    
    # mod.plot_QC()
    return adata_ref

def cl2_model(adata_ref,
            adata,
            N_cells_per_location=4,
            model_kwargs={}, 
            lasso_data=None, 
            file_name='deconv',
            accelerator="gpu",
            sp_model=None,
            verbose=True,
            save_path=None
            ):
    
    import cell2location
    import scvi
    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    # find mitochondria-encoded (MT) genes
    adata.var['mt_gene'] = [gene.startswith('mt-') for gene in adata.var_names]
    # remove MT genes for spatial mapping (keeping their counts in the object)
    adata = adata[:, ~adata.var_names.str.startswith('mt-')]

    if lasso_data is not None:
        adata = adata[lasso_data, :]

    intersect = np.intersect1d(adata.var_names, inf_aver.index)
    adata = adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata)

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=N_cells_per_location,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection (using default here):
        detection_alpha=20,
        **model_kwargs
        )
    
    mod.view_anndata_setup()
    mod.train(max_epochs=8000,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            accelerator="gpu")

    adata = mod.export_posterior(
        adata, sample_kwargs={
            'num_samples': 500, 
            'batch_size': mod.adata.n_obs, 
            'accelerator': 'gpu'}
    )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    # mod.plot_history(1000)
    # plt.legend(labels=['full data training']);

    if save_path is not None:
        mod.save(f"{save_path}", overwrite=True)
        adata.write(f'{file_name}.h5ad')
    
    return adata, mod

def run_cell2location(
    adata, 
    adata_ref, 
    batch_key='orig.ident',
    labels_key='subtype',
    sample_key=None,
    filtered=False,
    save_path='./', #TODO None
    sample='deconv', 
    lasso_data=None,
    sc_model=None, 
    sp_model=None,
    N_cells_per_location=2, 
    model_kwargs={},
    accelerator="gpu",
    verbose=True
    ):
    
    import scvi
    
    # filter the object
    if not filtered:
        if verbose: 
            print('filter single cell data...')
            
        from cell2location.utils.filtering import filter_genes
        selected = filter_genes(adata_ref, 
                                cell_count_cutoff=5,
                                cell_percentage_cutoff2=0.03,
                                nonz_mean_cutoff=1.12)
        adata_ref = adata_ref[:, selected].copy()
    
    if verbose: 
        print('Estimation of reference cell type signatures...')
    
    adata_ref = regression_model(adata_ref, 
                                save_path=save_path,
                                batch_key=batch_key,
                                labels_key=labels_key,
                                accelerator="gpu",
                                sc_model=sc_model
                                ) 
    
    if verbose: 
        print('training cell2location model...')
        
    adata, mod = cl2_model(adata_ref,
                    adata,
                    N_cells_per_location=N_cells_per_location, 
                    model_kwargs=model_kwargs,
                    lasso_data=lasso_data,
                    sp_model=sp_model
                    )
    
    return adata, mod

import sys
import argparse

parser = argparse.ArgumentParser(description="running cell2location")
parser.add_argument("-sp", metavar="file", required=True, help="spatial adata")
parser.add_argument("-ref", metavar="file", required=True, help="refdata of sc")
parser.add_argument("-outpfx", metavar="str", required=True, help="outpfx")
parser.add_argument("-label", metavar="str", required=True, help="label key")
parser.add_argument("-num", type=int, required=True, help="N_cells")
args = parser.parse_args()

adata_spatial = sc.read(args.sp)
adata_ref = sc.read(args.ref)
adata_ref.__dict__['_raw'].__dict__['_var'] = adata_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

adata_spatial.X =  adata_spatial.raw.X.copy()
adata_spatial.raw = adata_spatial

outpfx = args.outpfx
label_key = args.label
if adata_ref.raw is not None:
    adata_ref = adata_ref.raw.to_adata()

if "features" in adata_ref.var.keys():
    adata_ref.var.index = adata_ref.var["features"]
adata_ref.var["features"] = adata_ref.var.index
adata_spatial.obs["sample"] = outpfx

outdir = os.path.abspath(os.getcwd())
outdir = f"{outdir}/{outpfx}.{label_key}"
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.chdir(outdir)
adata_ref.obs[label_key] = adata_ref.obs[label_key].astype("category")
if "batch" not in adata_ref.obs.columns:
    adata_ref.obs["batch"] = adata_ref.obs["orig.ident"]

adata, mod = run_cell2location(adata_spatial, adata_ref, accelerator='gpu', labels_key=label_key, N_cells_per_location=args.num)

adata.obsm['stds_cell_abundance_w_sf'].to_csv(f'{outpfx}.denovo.stds_cell_abundance_w_sf.csv')
adata.obsm['means_cell_abundance_w_sf'].to_csv(f'{outpfx}.denovo.means_cell_abundance_w_sf.csv')
adata.obsm['q05_cell_abundance_w_sf'].to_csv(f'{outpfx}.denovo.q05_cell_abundance_w_sf.csv')
adata.obsm['q95_cell_abundance_w_sf'].to_csv(f'{outpfx}.denovo.q95_cell_abundance_w_sf.csv')
adata.write_h5ad(f'{outpfx}.denovo.h5ad')

###neighborhood 分析
# compute KNN using the cell2location output stored in adata.obsm
sc.pp.neighbors(adata, use_rep='q05_cell_abundance_w_sf',
                n_neighbors = 15)

# Cluster spots into regions using scanpy
sc.tl.leiden(adata, resolution=1.1)

# add region as categorical variable
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
# compute UMAP using KNN graph based on the cell2location output
sc.tl.umap(adata, min_dist = 0.3, spread = 1)
adata.write_h5ad(f'{outpfx}.denovo.h5ad')


###估计空间数据中每个基因的细胞类型特异性表达
# Compute expected expression per cell type
expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

# Add to anndata layers
for i, n in enumerate(mod.factor_names_):
    adata.layers[n] = expected_dict['mu'][i]

# Save anndata object with results
adata_file = f"{outpfx}.sp.denovo.expected_per_cell_type.h5ad"
adata.write(adata_file)
