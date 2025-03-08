recluster subcluster


usage: 20250308_recluster.r [-h] --input-rds INPUT_RDS  
                            [--output-dir OUTPUT_DIR] [--name NAME]  
                            [--cluster-id CLUSTER_ID] [--batch BATCH]  
                            [--method {integrate,harmony,none}] [--dims DIMS]  
                            [--resolution RESOLUTION]  

Process single-cell transcriptomics data with batch correction, clustering,
and differential expression analysis.

options:
  -h, --help            show this help message and exit
  --input-rds INPUT_RDS
                        Path to input Seurat RDS file (required)
  --output-dir OUTPUT_DIR
                        Output directory [default: current directory './']
  --name NAME           Prefix for output files [default:
                        'single_cell_analysis']
  --cluster-id CLUSTER_ID
                        Comma-separated cluster IDs to subset (e.g., '1,2,3'),
                        optional
  --batch BATCH         Batch variable for correction [default: 'replicate']
  --method {integrate,harmony,none}
                        Batch correction method: 'integrate' (CCA), 'harmony',
                        or 'none' [default: 'integrate']
  --dims DIMS           PCA dimensions for UMAP/clustering, R expression
                        (e.g., '1:20') [default: '1:15']
  --resolution RESOLUTION
                        Clustering resolution, numeric value [default: 0.3]
