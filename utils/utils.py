import pandas as pd
import numpy as np
import os

import utils_genelist as genelist

# Full donor ID mapping
full_donor_id_dict = {
    'HD4246': 'HD42466',
    'HD4249': 'HD42499',
    'HD4181': 'HD41818',
    'HD4241': 'HD42410',
    'HD4254': 'HD42541',
    'HD4309': 'HD43094',
    'HD4362': 'HD43627',
    'HD4236': 'HD42363',
    'HD4272': 'HD42722',
    'HD4277': 'HD42775',
}


# Set color code in anndata
def set_color_code(adata, color_code, feature=None):

    # Define for which features to set color
    if isinstance(feature, str):
        features = [feature]
    elif isinstance(feature, list):
        features = feature
    else:
        features = list(color_code.keys())

    # Define colors for features
    for feature in features:
        if feature not in adata.obs.columns:
            raise ValueError(f"Feature {feature} not found in adata.obs.")
        # Set category order in obs
        if not pd.api.types.is_categorical_dtype(adata.obs[feature]):
            adata.obs[feature] = adata.obs[feature].astype('category')
        # Order categories according to color code: would also need subsetting of color in uns slot etc
        # categs_avail = [cat for cat in color_code[feature]['order'] if cat in adata.obs[feature].cat.categories]
        adata.obs[feature] = adata.obs[feature].cat.reorder_categories(
            color_code[feature]['order'], ordered=True
        )
        # Set colors in uns
        key = f"{feature}_colors"
        adata.uns[key] = color_code[feature]['colors']

    return adata


def set_uns_colors(adata, color_code, key=None):
    """
    Sets adata.uns[f"{key}_colors"] to a list of hex colors based on color_code,
    compatible with scanpy plotting functions.
    
    Parameters:
    - adata: AnnData object
    - key: str, name of obs column
    - color_code: dict with keys 'order' and 'colors', e.g.:
        {
            'order': ['cat1', 'cat2'],
            'colors': {
                'cat1': '#hex1',
                'cat2': '#hex2'
            }
        }
    """

    # Handle single key or nested dictionary
    if key is None:
        keys = color_code.keys()
    else:
        keys = [key]

    # Iterate over color keys
    for key in keys:

        if key not in adata.obs.columns:
            print(f"[WARNING] Color key '{key}' not in adata.obs.columns. Skipping...")
            continue

        # Extract order and color mapping
        order = color_code[key]['order']
        colors = color_code[key]['colors']

        # Subset to categories present in adata.obs[key]
        # Set category order in obs
        if not pd.api.types.is_categorical_dtype(adata.obs[key]):
            adata.obs[key] = adata.obs[key].astype('category')
        present_cats = [cat for cat in order if cat in adata.obs[key].cat.categories]
        if len(present_cats) < len(order):
            order = present_cats
            colors = {cat: colors[cat] for cat in present_cats}

        # Ensure obs[key] is categorical with the correct order
        adata.obs[key] = pd.Categorical(adata.obs[key], categories=order, ordered=True)

        # Generate color list in correct order
        try:
            color_list = [colors[cat] for cat in order]
        except KeyError as e:
            raise KeyError(f"Missing color for category '{e.args[0]}' in '{key}' color definition.")

        # Assign to .uns
        adata.uns[f"{key}_colors"] = color_list

    # return adata


### SC status
def add_SC_status(adata, species = "Hs"):
    # Define gene sets and colors
    if species == "Hs":
        features = ["LGR5", "TACSTD2"]
    else:
        features = ["Lgr5", "Tacstd2"]
    colors = ['green', "red"]
    # For each cell check if feature gene expression is > 0
    for feature, color in zip(features, colors):
        if feature in adata.var_names:
            adata.obs[f"{feature}_pos"] = adata.obs_vector(feature, layer = "counts") > 0
        else:
            print(f"Feature {feature} not found in dataset.")
    ### Assign SC population status using gene expression levels
    conditions = [
        (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == True),
        (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == False),
        (adata.obs[f"{features[0]}_pos"] == False) & (adata.obs[f"{features[1]}_pos"] == True),
    ]
    # Corresponding choices for each condition
    choices = ["both", features[0], features[1]]
    # Create the new column based on the conditions and choices
    adata.obs['SC_status'] = np.select(conditions, choices, default="none")
    return adata


### EMP1 / TROP2 status
def add_EMP1_status(adata, species = "Hs"):
    # Define gene sets and colors
    if species == "Hs":
        features = ["EMP1", "TACSTD2"]
        colors = ['blue', "green"]
    else:
        features = ["Emp1", "Tacstd2"]
        colors = ['blue', "green"]
    # For each cell check if feature gene expression is > 0
    for feature, color in zip(features, colors):
        if feature in adata.var_names:
            adata.obs[f"{feature}_pos"] = adata.obs_vector(feature, layer = "counts") > 0
        else:
            print(f"Feature {feature} not found in dataset.")

    ### Assign population status using gene expression levels
    conditions = [
        (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == True),
        (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == False),
        (adata.obs[f"{features[0]}_pos"] == False) & (adata.obs[f"{features[1]}_pos"] == True),
    ]
    # Corresponding choices for each condition
    choices = ["both", features[0], features[1]]
    # Create the new column based on the conditions and choices
    adata.obs['EMP1_status'] = np.select(conditions, choices, default="none")
    return adata



### Gene set related functions

# Get gene sets
def get_gene_sets(adata, verbose=True):

    # Initiate gene set collection
    gene_sets = []

    ### YAP score (Gregorieff, Nature, 2015)
    if verbose:
        print("Get Gregorieff YAP gene set")
    score_name = "score_YAP"
    gene_list = genelist.get_Gregorieff_2015(species="Hs", genes_only = True, verbose=verbose)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list, verbose=verbose)

    ### Fetal score (Mustata, CellRep, 2013)
    if verbose:
        print("Process fetal gene set")
    score_name = "score_fetal"
    gene_list = genelist.get_Mustata_2013(species="Hs", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list, verbose=verbose)

    ### Lgr5 score (Muñoz, EMBOJ, 2011)
    if verbose:
        print("Process LGR5 gene set")
    score_name = "score_Lgr5"
    gene_list = genelist.get_Munoz_2011(species="Hs", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list, verbose=verbose)

    ### Ganesh-related gene sets
    if verbose:
        print("Process Ganesh gene set collection")
    gene_list = genelist.get_Moorman_2023(extract_geneset="all", genes_only = True)
    for score_name in gene_list['geneset'].unique():
        df = gene_list[gene_list['geneset'] == score_name]
        genes = df['gene'].to_list()
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes, verbose=verbose)
        
    ### Emp1 score (Cañellas-Socias, Nature, 2022)
    if verbose:
        print("Process EMP1-related gene sets")
    scores = {
        "score_allHR": 'All-HR',
        "score_coreHR": 'coreHRC',
        "score_epiHR": 'EpiHR',
        "score_tmeHR": 'TME-HR',
        "score_Lgr5_Batlle": 'Lgr5 signature',	
        "score_Wnt_Batlle": 'Wnt ON signature (Morral et al. 2020)',
        "score_mKi67_Batlle": 'mKi67 (Basak et al. 2014)',
        "score_Yap_Batlle": 'YAP_22 (Wang et al. 2018)'
    }
    gene_list = genelist.get_CanellasSocias_2022(species="Hs")
    # Iterate over all scores
    for score_name,col_name in scores.items():
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list[col_name], verbose=verbose)

    ### Batlle 2025 paper (Centonze et al., Cancer Discovery, 2025) -- added 01.05.2026
    if verbose:
        print("Process Centonze (Batlle) 2025 gene sets")
    scores = {
        "score_Revival_SCs_(Vazquez)": "Revival Stem Cells (Vazquez et al. 2022)",
        "score_Revival_SCs_(Ayyaz)": "Revival Stem Cells (SSC2c) (Ayyaz et al. 2019)",
        "score_Immature_enterocytes_(Smillie)": "Immature enterocytes (Smillie et al. 2019)",
        "score_Basal_PDAC_(Raghavan)": "Basal pancreatic cancer cells (Raghavan et al. 2021)",
        "score_Proliferation_(Merlos)": "Proliferation (Merlos et al. 2011)",
        # "score_iCMS2_(Joanito)": "iCMS2 (Joanito et al. 2022)",
        # "score_iCMS3_(Joanito)": "iCMS3 (Joanito et al. 2022)",
    }
    gene_list = genelist.get_Centonze_2025()
    ### Convert to desired species
    # data frame to long table format
    gene_list_long = gene_list.melt(var_name='geneset', value_name='gene').dropna()
    # Remove unwanted spaces from gene names
    gene_list_long['gene'] = (
        gene_list_long['gene']
        .astype(str)
        .str.strip()  # remove normal spaces, tabs, newlines
        # .str.replace(r'[\u00A0\u200B\u200C\u200D\uFEFF]', '', regex=True)  # remove non-breaking / zero-width spaces
    )    
    # Detect gene format per geneset
    geneset_species = (
        gene_list_long.groupby('geneset')['gene']
        .apply(lambda x: 'Hs' if (x.str.isupper().mean() > 0.75) else 'Mm')
    )
    # Split gene sets based on detected origin
    human_sets = geneset_species[geneset_species == 'Hs'].index
    mouse_sets = geneset_species[geneset_species == 'Mm'].index
    gene_list_human = gene_list_long[gene_list_long['geneset'].isin(human_sets)]
    gene_list_mouse = gene_list_long[gene_list_long['geneset'].isin(mouse_sets)]
    # print(f"Detected {len(human_sets)} human-like gene sets and {len(mouse_sets)} mouse-like gene sets")
    # Perform gene conversion
    if verbose:
        print(f"Converting {len(gene_list_mouse)} mouse genes to human via capitalization")
    gene_list_converted = gene_list_mouse.copy()
    gene_list_converted['gene'] = gene_list_converted['gene'].str.upper()
    gene_list_long = pd.concat([gene_list_human, gene_list_converted], ignore_index=True)
        
    ### Iterate over all scores
    for score_name, col_name in scores.items():
        genes = gene_list_long[gene_list_long['geneset'] == col_name]['gene'].to_list()
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes, verbose=verbose)


    ### Batlle-related gene set collection of Nuria V
    if verbose:
        print("Process Batlle-related gene set collection of Nuria V")
    gene_list = genelist.get_collection_NV()
    # Iterate over all gene sets
    for score_name in gene_list.columns:
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, 'score_'+score_name, gene_list[score_name], verbose=verbose)    

    ### Flatten to one gene / row dataframe
    gene_sets = pd.DataFrame([
        {**d, 'gene': gene} for d in gene_sets for gene in d['genes']
    ]).drop('genes', axis=1)
    # Reorder columns
    gene_sets = gene_sets[['gene_set', 'gene', 'weight']]

    # Add results to anndata slot
    adata.uns['gene_sets'] = gene_sets

    return adata



### Import obs output
def import_obs_features(adata, analysis, file_id="", filename=None, path=None, if_avail=False, dtype={}):
    
    # Define path of obs file
    if path is None:
        path = os.path.join("../data/", analysis)
    
    ### Read in obs file
    file = filename if filename else f"{file_id}_obs_{analysis}.csv"
    # Check if file exists
    if not os.path.exists(os.path.join(path, file)):
        if if_avail:
            print(f"Obs file for {analysis} not found. Skipping import...")
            print(os.path.join(path, file))
            return adata
        else:
            raise FileNotFoundError(f"File {file} not found in {path}.")
    df = pd.read_csv(os.path.join(path, file), index_col=0,dtype=dtype)

    # Merge obs file with adata
    df_col_new = [col not in adata.obs.columns for col in df.columns]
    if np.all(df_col_new):
        print(f"{sum(df_col_new)} new obs features will be added.")
    elif np.any(df_col_new):
        print(f"Following columns are already part of adata.obs: {df.columns[~pd.Series(df_col_new)]}.")
        print(f"{sum(df_col_new)} new obs features will be added.")
        df = df.loc[:, df_col_new]
    else:
        print("All columns are already part of adata.obs.")
        return adata

    adata.obs = adata.obs.merge(df, left_index=True, right_index=True, how='left')
    return adata


def import_uns_feature(adata, analysis, file_id="", filename=None, path=None, if_avail=False):
    
    # Define path of obs file
    if path is None:
        path = os.path.join("../data/", analysis)
    
    # Read in obs file
    file = filename if filename else f"{file_id}_gene_set_collection.csv" if analysis == 'gene_sets' else f"{file_id}_uns_{analysis}.csv"
    # Check if file exists
    if not os.path.exists(os.path.join(path, file)):
        if if_avail:
            print(f"Uns file for {analysis} not found. Skipping import...")
            return adata
        else:
            raise FileNotFoundError(f"File {file} not found in {path}.")
    df = pd.read_csv(os.path.join(path, file), index_col=None)

    # Add uns df to adata
    if analysis not in adata.uns:
        adata.uns[analysis] = df
    else:
        # Print warning if uns will be replaced
        if not adata.uns[analysis].equals(df):
            print(f"Warning: {analysis} already part of adata.uns. Will be replaced.")
            adata.uns[analysis] = df
        else:
            print(f"Analysis {analysis} already part of adata.uns. Skipping...")
    
    return adata


def import_magic(
    adata,
    final_id=None,
    path_magic="../data/magic/",
    filename_magic=None,
    layer_key="magic",
    fill_missing_genes=False,
    verbose=True,
):
    """
    Import MAGIC-imputed expression matrix from CSV and store in adata.layers[layer_key].

    Designed for AnnData layers:
      - Output MUST match adata.shape (n_obs, n_vars).
      - MAGIC matrix may come from a larger dataset:
          * MAGIC can have extra cells/genes (they will be dropped, reported).
          * But MAGIC must contain ALL adata.obs_names.
      - Missing genes can optionally be added as zeros (fill_missing_genes=True).

    Parameters
    ----------
    adata : AnnData
        Target AnnData object (not modified except adding/updating the layer).
    final_id : str, optional
        Used to derive filename if filename_magic is not provided.
    path_magic : str
        Folder containing MAGIC CSV.
    filename_magic : str, optional
        Explicit filename of MAGIC CSV. If None, derived from final_id.
    layer_key : str
        Name of the layer to store MAGIC values in.
    fill_missing_genes : bool
        If True, genes present in adata.var_names but missing from MAGIC are added with zeros.
        If False, missing genes raise an error.
    verbose : bool
        Print informational messages.

    Returns
    -------
    adata : AnnData
        Same object, with adata.layers[layer_key] set.
    """

    # Derive filename
    if final_id is None and filename_magic is None:
        raise ValueError("Either final_id or filename_magic must be provided.")
    if filename_magic is None:
        filename_magic = f"{final_id}_magic.csv"

    fpath = os.path.join(path_magic, filename_magic)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"MAGIC file not found: {fpath}")

    if verbose:
        print(f"Importing MAGIC results from: {fpath}")

    magic = pd.read_csv(fpath, index_col=0)

    # --- Enforce the requested rule: MAGIC dataset must be >= adata in obs dimension
    if magic.shape[0] < adata.n_obs:
        raise ValueError(
            f"MAGIC has fewer observations than adata "
            f"({magic.shape[0]} < {adata.n_obs}). "
            f"This importer only supports importing into adata subsets of MAGIC (MAGIC n_obs >= adata n_obs)."
        )

    # --- Check that all adata obs exist in MAGIC (subset import)
    missing_obs = adata.obs_names.difference(magic.index)
    if len(missing_obs) > 0:
        raise ValueError(
            f"{len(missing_obs)} adata observation(s) are missing from MAGIC matrix. "
            f"Cannot create a layer without values for every adata.obs_names entry. "
            f"Example missing obs: {list(missing_obs[:5])}"
        )

    # --- Report extra obs that will be dropped
    extra_obs = magic.index.difference(adata.obs_names)
    if verbose and len(extra_obs) > 0:
        print(f"Dropping {len(extra_obs)} extra observation(s) from MAGIC (keeping adata subset).")

    # --- Handle genes
    missing_genes = adata.var_names.difference(magic.columns)
    if len(missing_genes) > 0:
        if fill_missing_genes:
            if verbose:
                print(
                    f"{len(missing_genes)} gene(s) in adata.var_names are missing from MAGIC. "
                    f"Adding them as zeros (fill_missing_genes=True)."
                )
            # Add missing gene columns as zeros
            magic = magic.reindex(columns=magic.columns.union(adata.var_names), fill_value=0)
        else:
            raise ValueError(
                f"{len(missing_genes)} gene(s) in adata.var_names are missing from MAGIC. "
                f"Set fill_missing_genes=True to pad them with zeros. "
                f"Example missing genes: {list(missing_genes[:10])}"
            )

    # --- Report extra genes that will be dropped
    extra_genes = magic.columns.difference(adata.var_names)
    if verbose and len(extra_genes) > 0:
        print(f"Dropping {len(extra_genes)} extra gene(s) from MAGIC (keeping adata subset genes).")

    # --- Subset + reorder EXACTLY to adata (layer-safe)
    magic = magic.loc[adata.obs_names, adata.var_names]

    # --- Final safety checks (should always pass now)
    if magic.shape != adata.shape:
        raise ValueError(f"Final MAGIC matrix shape {magic.shape} != adata.shape {adata.shape} (unexpected).")
    if not np.array_equal(magic.index, adata.obs_names):
        raise ValueError("Final MAGIC obs order does not match adata.obs_names (unexpected).")
    if not np.array_equal(magic.columns, adata.var_names):
        raise ValueError("Final MAGIC var order does not match adata.var_names (unexpected).")

    # --- Assign layer
    adata.layers[layer_key] = magic.to_numpy()
    
    if verbose:
        print(f"Stored MAGIC matrix in adata.layers['{layer_key}'] with shape {adata.layers[layer_key].shape}")
        
    # Free up memory
    del magic
