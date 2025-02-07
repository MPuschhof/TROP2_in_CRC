import os, re
import pickle
import pandas as pd
import scanpy as sc
from pathlib import Path

    
def set_paths(DATA_PATH=None, OUT_DATA_PATH=None, PLOT_PATH=None, name=None,
              subdir_outdata=True, subdir_plot=True):
    """
    This function sets the paths
    """
    S_PATH = "/".join(os.path.realpath(__file__).split(os.sep)[:-1])
    
    ### Adapt paths if provided through arg parse
    # Input data path
    DATA_PATH = os.path.join(S_PATH, "../data") if DATA_PATH is None else DATA_PATH
    # Output data path
    if subdir_outdata:
        OUT_DATA_PATH = os.path.join(S_PATH, "../data", name) if OUT_DATA_PATH is None else DATA_PATH
        Path(OUT_DATA_PATH).mkdir(parents=True, exist_ok=True)
    else:
        OUT_DATA_PATH = os.path.join(S_PATH, "../data") if OUT_DATA_PATH is None else OUT_DATA_PATH
    # Output plot path
    if subdir_plot:
        PLOT_PATH =  os.path.join(S_PATH, "../plots", name) if PLOT_PATH is None else PLOT_PATH
        Path(PLOT_PATH).mkdir(parents=True, exist_ok=True)
    else:
        PLOT_PATH = os.path.join(S_PATH, "../plots") if PLOT_PATH is None else PLOT_PATH
    
    return S_PATH, DATA_PATH, OUT_DATA_PATH, PLOT_PATH
    

# def save_to_subdir(parent_path, subdir=None):
#     if subdir is None:
#         warnings.warn("No name for subdirectory provided. Will be saved in main directory.")
#         return
#     else:
#         subdir_path = Path(parent_path, subdir)
#         subdir_path.mkdir(parents=True, exist_ok=True)
#         return subdir_path

# def read_h5_from_crg_multi(dir_data, sample_name):
#     adata = sc.read_10x_h5(os.path.join(dir_data, sample_name, 'count/sample_filtered_feature_bc_matrix.h5'))
#     adata.var_names_make_unique()
#     return adata

# def read_raw_sample_barcodes(dir_data, sample_name):
#     bc_data = pd.read_csv(os.path.join(dir_data, sample_name, 'count/sample_filtered_barcodes.csv'), header = None)
#     bc_data.rename(columns={0: 'ref_genome', 1: 'cell_bc'}, inplace=True)
#     return bc_data

# def add_metadata_to_anndata(adata, meta_data):
#     for col in meta_data.columns:
#         adata.obs[col] = pd.Series(meta_data[col], dtype="string").values
#     # ensure all obs columns are of dtype 'str' as otherwise saving to h5ad might fail
#     for col in adata.obs.columns:
#         adata.obs[col] = adata.obs[col].astype(str)
#     return adata

# def lowercase_mt_metric_cols(adata):
#     # Rename columns in `obs` that end with "_MT" to end with "_mt"
#     adata.obs.rename(columns={col: col.replace("_MT", "_mt") for col in adata.obs.columns if col.endswith("_MT")}, inplace=True)

#     # Check for a column named exactly "MT" in `var` and rename it to "mt"
#     if "MT" in adata.var.columns:
#         adata.var.rename(columns={"MT": "mt"}, inplace=True)
#     return adata

# def read_h5ad(subdir, identifier):
    
#     # Identify correct file
#     dir_data = os.path.join("../data/", subdir)
#     files = os.listdir(dir_data)
    
#     reg = re.compile(rf'.*{re.escape(identifier)}.*\.h5ad') # Compile the regex
#     data = list(filter(reg.search, files)) 
#     adata = sc.read_h5ad(os.path.join(dir_data, data[0]))

#     return adata

# ### More MP functions
def find_key_for_value(my_dict, search_value): # for unique results
    flag = True
    for key in my_dict.keys():
        if search_value in my_dict[key]:
            flag = False
            return key
    if flag:
        print(f"Warning: {search_value} not found in dictionary")
        return False


# ### Pseudotime-related functions
# ### RELATED TO SCVELO
# def reformat_loom_file(dir_loom, file):
    
#     # Try to open loom file
#     try:
#         with loompy.connect(os.path.join(dir_loom, file), validate=True) as ds:
#             pass
#     except ValueError:
#         #print("File format probably not fulfilling to loompy / scvelo expectations.")
#         reformat_loom = True
#     else:
#         print("File can be opened.")
#         return
        
#     # Adapt format to the one expected by loompy and scvelo
#     if reformat_loom:
#         with loompy.connect(os.path.join(dir_loom, file), validate=False) as ds:
#             for attr in ['Accession', 'Chromosome', 'Gene', 'Strand', 'CellID']:
#                 if attr in ds.ra:
#                     #print(f"{attr} found in ds.ra, convert to str")
#                     ds.ra[attr] = ds.ra[attr].astype(str)
#                 if attr in ds.ca:
#                     #print(f"{attr} found in ds.ca, convert to str")
#                     ds.ca[attr] =  ds.ca[attr].astype(str)
                    
#         # Try to open loom file again
#         try:
#             with loompy.connect(os.path.join(dir_loom, file), validate=True) as ds:
#                 pass
#         except ValueError:
#             print("File format still not fulfilling to loompy / scvelo expectations.")
#             return
#         else:
#             print("File successfully reformatted.")
#             return
    
# def harmonize_bc_names(ldata, sample_ID, bc_dict):
#     # Test if ":" is still part of index column as in original loom file    
#     try:
#         bc_ldata = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
#     except IndexError:
#         print("Barcode names might already be harmonized. Please check")
#         return
#     else:
#         # Modify barcode names
#         bc_ldata = [bc[0:len(bc)-1] + "-1" for bc in bc_ldata]
#         ldata.obs.index = bc_ldata
         
