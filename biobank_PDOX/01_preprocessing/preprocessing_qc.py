import os
import scanpy as sc
import numpy as np
import scanpy.external as sce
import argparse

import utils, utils_spec




#################### Argument parsing and related parameters

# # Define input parameters
# DATA_PATH = '../data/crg_multi/'
# analysis = 'pp_qc'
# sample_id = ''

# # Derive output paths
# OUT_DATA_PATH = os.path.join(DATA_PATH, f'../{analysis}')
# PLOT_PATH = os.path.join(DATA_PATH, f'../../plots/{analysis}')

# # Generate output paths if non-existant
# Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
# Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)

# Define argument parser
parser = argparse.ArgumentParser(prog='qc', description='Sample-wise QC and filtering')
parser.add_argument('-i', '--input_path', help='path to the input data', required=True)
parser.add_argument('-f', '--filename', help='filename of input data', required=True)
parser.add_argument('-an', '--analysis_name', help='analysis name', required=True)
parser.add_argument('-o', '--output_path', help='path for data output', required=False)
parser.add_argument('-p', '--plot_path', help='path for figure output', required=False)

# Extract values
args = vars(parser.parse_args())
dir_data = args['input_path']
filename = args['filename']
analysis = args['analysis_name']
dir_out = args['output_path']
dir_plot = args['plot_path']

# Get necesary paths and create folders if necessary   
S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH = utils.set_paths(dir_data, dir_out, dir_plot, analysis)


#################### Input parameters for processing

from utils_spec import mt_thr, gene_qnt, min_genes, doublet_thr, min_cells

#################### Derive further  variables from input

sample_id = filename.split("_filtered_feature_bc_matrix")[0]


#################### Functions and code to execute
        
def preprocess_and_qc(adata):

    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.to_series().str.lower().str.contains("mt-")
    if adata.var['mt'].sum() == 0:
        print("Warning: no mitochondrial genes detected. Check gene nomenclature.")
        
    # QC metric: gene counts
    sc.pp.calculate_qc_metrics(adata, expr_type = 'counts', 
                               var_type = 'genes', qc_vars=["mt"], inplace=True)
    
    # calculate gene threshold, distribution plot generated later
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt)

    # QC metric: doublets
    sce.pp.scrublet(adata, verbose = False) # warning on making a copy to be expected
        
    # Filtering data
    print(f"AnnData shape before any filtering: {np.shape(adata.X)}")
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes) # defined in plotting section
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells) # lower limit
    gene_thr = np.quantile(adata.obs.n_genes_by_counts, gene_qnt) # upper limit
    adata = adata[adata.obs.n_genes_by_counts < gene_thr, : ]
    # Filter mitochondrial content
    adata = adata[adata.obs.pct_counts_mt < mt_thr] # defined in plotting section
    
    # Filter doublets
    adata = adata[adata.obs.doublet_score < doublet_thr, : ]
    print(f"AnnData shape after filtering: {np.shape(adata.X)}")
        
    # Save processed file to disk
    print("Write anndata file to disk")
    adata.write(os.path.join(OUT_DATA_PATH, f"{sample_id}_filtered.h5ad"))
        

                
        
if __name__ == "__main__":
    
    print(f"Process sample {filename}")
    
    # Read h5 files from cellranger-multi output
    adata = sc.read_10x_h5(os.path.join(DATA_PATH, filename))
    adata.var_names_make_unique()

    # Add experimental meta data
    adata = utils_spec.get_meta_data_C01a006(adata, sample_id)
        
    # Run scRNA-seq QC with given parameters
    preprocess_and_qc(adata)