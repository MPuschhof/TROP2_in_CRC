import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import matplotlib.pyplot as plt
import scanpy.external as sce
import seaborn as sns
import re
import math
import warnings

dir_gene_lists = '/home/m014f/gene_lists/gene_lists_CRC'


##### Generic functions for gene set handling

def remove_nan_cols(df):
    cols_to_remove = []
    for col in df.columns:
        if df[col].isna().all():
            cols_to_remove.append(col)
            df.drop(columns=[col], inplace=True)
    if cols_to_remove:
        warnings.warn(f"Removed columns containing only NaNs: {', '.join(cols_to_remove)}")
    return df

    
def filter_genes(geneset, FC='FC', thres=None, absolute = False):
    if thres is None:
        raise ValueError("Threshold value must be defined")
    if FC not in geneset.columns:
        raise AttributeError("Column to filter on not defined")
    
    # Return reduced dataframe
    if absolute:
        return geneset[abs(geneset[FC]) > thres]
    else:
        return geneset[geneset[FC] > thres]
    

def add_to_gene_sets(adata, gene_sets, score_name, genes, weight=1, verbose=True):
    # Called in utils_spec.get_gene_sets()
    
    # Remove NaNs and duplicates
    genes = [x for x in genes if not (isinstance(x, float) and math.isnan(x))]
    genes = list(set(genes))

    # Subset to genes detected in this dataset
    N_geneset = len(genes)
    N_incl = pd.Series(genes).isin(adata.var.index).sum()
    if N_incl > 0:
        if verbose:
            print(f'{N_incl/N_geneset*100:.2f}% of the geneset {score_name} ({N_incl} out of {N_geneset} genes) are included in dataset.')
        # Create a dictionary for the current gene set and append it to the list
        temp_gene_set = {
            "gene_set": score_name,
            "genes": [gene for gene in genes if gene in adata.var.index],
            "weight": weight
        }
        gene_sets.append(temp_gene_set)
    else:
        if verbose:
            print(f'{N_incl/N_geneset*100:.2f}% of the geneset {score_name} ({N_incl} out of {N_geneset} genes) are included in dataset.')
        raise UserWarning(f"Gene set {score_name} will be skipped.") 
    
    return gene_sets


def compute_gene_score(adata, gene_list, score_name, use_raw=True, layer = None, verbose=True): 
    # Note: changed 11.05.2026 to not copy unnecessary counts.
    # Implementation: do not use return value.

    if not isinstance(adata, anndata.AnnData):
        raise TypeError("adata must be of type AnnData")
    
    # Remove NaN from gene list if necessary
    gene_list = [x for x in gene_list if not (isinstance(x, float) and math.isnan(x))]
    
    # Reduce list to non-duplicated entries if necessary
    gene_list = list(set(gene_list))
    
    if layer:
        if layer not in adata.layers:
            raise AttributeError(f"Layer {layer} not found in adata.")

        use_raw=False
        if verbose:
            print(f"Running score_genes on layer {layer} (use_raw set to False).")
        # X_orig = adata.X.copy()
        # adata.X = adata.layers[layer].copy()

    elif use_raw and adata.raw is not None:
        if verbose:
            print("Using adata.raw for gene score computation.")

    else:
        raise AttributeError("Data layer not found in adata.")

    # Generate / Print gene stats
    N_geneset = len(gene_list)
    N_incl = pd.Series(gene_list).isin(adata.var.index).sum()
    if N_incl > 0:
        #print(f'{N_incl/N_geneset*100:.2f}% of the geneset {score_name} ({N_incl} out of {N_geneset} genes) are included in dataset.')
        if verbose:
            print(f'Computing gene score for {score_name}')
        adata = sc.tl.score_genes(
            adata, 
            gene_list=gene_list, 
            score_name=score_name, 
            copy=False, 
            use_raw=use_raw,
            layer=layer,
        )
    else:
        if verbose:
            print(f'{N_incl/N_geneset*100:.2f}% of the geneset {score_name} ({N_incl} out of {N_geneset} genes) are included in dataset.')
        raise UserWarning(f"Gene set {score_name} cannot be computed.") 
    
    # # Swap active layers
    # if not use_raw:
    #     if verbose:
    #         print("Restoring original matrix X.")
    #     adata.X = X_orig

    # return adata


def strip_gene_list(gene_list):
    '''
    # Remove NaN from gene list
    gene_list = gene_list[gene_list.notnull()]
    
    # Return stripped gene names
    return [gene.strip() for gene in gene_list]
    '''

    # Remove NaN from gene list
    gene_list = gene_list[gene_list.notnull()]
    
    # Initialize new list
    genes_stripped = []
    
    for gene in gene_list:
        # Check if gene is a string
        if isinstance(gene, str):
            genes_stripped.append(gene.strip())
        else:
            # For non-string values, you can choose to append them as-is,
            # convert them to a string, or handle them differently
            genes_stripped.append(gene)
    
    return genes_stripped

    

##### Functions related to individual gene sets

def get_Gregorieff_2015(species="Mm", genes_only=False, FC=1.2, verbose=True):
    
    # Gregorieff et al, Nature, 2015
    # Data as listed in ST2
    # Two YAP perturbation experiments performed in mice
    # (1) YAP KO vs control
    # (2) YAP overexpression (dox-inducible)

    # Valid input parameters
    if species not in ['Mm', 'Hs']:
        raise AttributeError("Valid species are {Mm, Hs}.")

    files = os.listdir(dir_gene_lists)

    # File identification
    if species == "Mm":
        reg = re.compile(rf'.*Gregorieff.*ST2\.csv')
    elif species == "Hs":
        reg = re.compile(rf'.*Gregorieff.*Hs\.csv')

    lists = list(filter(reg.search, files)) 
    geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]), skiprows=0)
    
    # Rename columns specific to mouse list import
    if species == "Mm":
        colnames = {
            'Gene_Symbol': 'gene',
            'GeneID': 'geneID',
        }
        geneset.rename(columns=colnames, inplace=True)
        
    # Rename columns specific to human list import
    elif species == "Hs":
        geneset.columns = geneset.columns.str.replace('.', ' ')
        colnames = {
            'YapTg Dox ': 'YapTg Dox-',
            'YapTg Dox  1': 'YapTg Dox+',
        }
        geneset.rename(columns=colnames, inplace=True)
        
    # Rename columns for all species
    colnames = {
        'Control': 'YapKO_wt',
        'Yap KO': 'YapKO_KO',
        'YapTg Dox-': 'YapOE_cntrl',
        'YapTg Dox+': 'YapOE_dox',
        'Combined fold change': 'fc_import',
        'Fold_log2': 'fc_log2_import',
    }
    geneset.rename(columns=colnames, inplace=True)

    # Compute fold change for each perturbation
    geneset['FC_YapKO'] = geneset['YapKO_KO']/geneset['YapKO_wt']
    geneset['FC_YapOE'] = geneset['YapOE_dox']/geneset['YapOE_cntrl']

    # Compute directionality of perturbation
    geneset['Corr_YapKO'] = geneset['FC_YapKO'] > 1
    geneset['Corr_YapOE'] = geneset['FC_YapOE'] > 1

    # Compute combined (absolute) fold change
    geneset['FC'] = (1/geneset['FC_YapKO'] + geneset['FC_YapOE'])/2
    geneset['FC_log2'] = np.log2(geneset['FC'])
    geneset['FC_log2_abs'] = abs(geneset['FC_log2'])
    
    if not genes_only:
        return geneset   

    else:
        # Filter out genes with opposing directionality
        if verbose:
            print(f"Genes with opposing directionality in YapKO and YapOE will be filtered out.")
        geneset['same_direction'] = (~ geneset['Corr_YapKO'] == geneset['Corr_YapOE'])
        geneset = geneset[geneset['same_direction']]
        
        # Filter on effect size
        if verbose:
            print(f"Genes with an absolute fold change of log2({FC}) will be kept.")
        geneset = filter_genes(geneset, FC = 'FC_log2_abs', thres = np.log2(FC))
        
        if species == "Mm":
            return geneset['gene'][geneset['Corr_YapOE']]
        elif species == "Hs":
            return geneset['Human_symbol'][geneset['Corr_YapOE']]
        
    
def get_Munoz_2011(species="Mm", genes_only=False):
    
    # Munoz, EMBO Journal, 2011
    # Data as listed in ST3 (gene names only)
    # Underlying comparison: Lgr5-high vs Lgr5-low cells
    
    # Valid input parameters
    if species not in ['Mm', 'Hs']:
        raise AttributeError("Valid species are {Mm, Hs}.")

    files = os.listdir(dir_gene_lists)

    # File identification and reading
    if species == "Mm":
        reg = re.compile(rf'.*Munoz.*ST3\.csv')
        lists = list(filter(reg.search, files)) 
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]), skiprows=2)
    elif species == "Hs":
        reg = re.compile(rf'.*Munoz.*Hs\.csv')
        lists = list(filter(reg.search, files)) 
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]))

    # Rename columns specific to mouse list import
    if species == "Mm":
        colnames = {
            'Official Gene Symbol': 'gene',
            'Overlap Group': 'enrichment_group',
            'Gene Symbol (Affymetrix)': 'gene_Affymetrix',
            'Gene Symbol (Agilent)': 'gene_Agilent',
        }
        # Only keep columns with non-NA values
        geneset = geneset.iloc[:, :4]
        
    # Rename columns specific to human list import
    elif species == "Hs":
        colnames = {
            'Official.Gene.Symbol': 'gene',
            'Overlap.Group': 'enrichment_group',
            'Gene.Symbol..Affymetrix.': 'gene_Affymetrix',
            'Gene.Symbol..Agilent.': 'gene_Agilent',
        }
    geneset.rename(columns=colnames, inplace=True)

    # Return gene list or gene table
    if genes_only:
        if species == "Mm":
            return geneset['gene']
        elif species == "Hs":
            return geneset['Human_symbol']
    else:
        return geneset


def get_Moorman_2023(extract_geneset="fetal", species="Hs", genes_only=False):
    
    # Moorman et al. Nature 2023 (@Ganesh)
    # Data as listed in ST3 (Hotspot factors) for squamous and neuroendocrine signatures
    # Data as listed in ST5 with log2FC for fetal sign
    
    # Valid input parameters
    if extract_geneset not in ['fetal', 'squamous', 'neuroendocrine', 'all']:
        raise AttributeError("Valid genesets are {fetal, squamous, neuroendocrine, all}.")
    if species not in ['Hs']:
        raise AttributeError("Only implemented for human genes.")

    # Gene set names
    gene_set_names = {
        'fetal': 'score_fetal_Moorman',
        'squamous': 'score_squamous_Moorman',
        'neuroendocrine': 'score_neuroendocrine_Moorman',
    }

    # Get list of all available files
    files = os.listdir(dir_gene_lists)
    files = [f for f in files if f.startswith('Moorman-Nature-2023') and f.endswith('.csv')]
    
    # Check for correct gene set file list
    if extract_geneset == 'all':
        if len(files) != 2:
            raise ValueError(f"{len(files)} files found for given geneset (expected 2).")
    else:
        # Filter for specific input file
        if extract_geneset == "fetal":
            files = [f for f in files if 'ST5' in f]
        else:
            files = [f for f in files if 'ST3' in f]
        if len(files) != 1:
            raise ValueError(f"{len(files)} files found for given geneset (expected 1).")

    # Define output column names
    colnames = {
        'Gene': 'gene',
        'Pval (FDR Adj)': 'p_val_FDR_adj',
        # 'R, KG146': 'Pearson_R_donor_KG146',
        # 'R, KG182': 'Pearson_R_donor_KG182',
        # 'R, KG150': 'Pearson_R_donor_KG150',
    }

    # Read human gene list and clean up
    genesets = []
    for file in files:
        df = pd.read_csv(os.path.join(dir_gene_lists, file))
        # Reduce to extract_geneset
        if 'ST3a' in file:
            keep_geneset = ['Squamous', 'Neuroendocrine'] if extract_geneset == 'all' else [extract_geneset.capitalize()]
            df = df[df['Annotation'].isin(keep_geneset)]
            df.reset_index(drop=True, inplace=True)
        # Rename columns
        df.rename(columns=colnames, inplace=True)
        # Remove rows with only NaN values
        df.dropna(how='all', inplace=True)
        # Add gene set name
        if (len(files) > 1):
            if file.endswith('fetal.csv'):
                df.insert(0, 'geneset', gene_set_names['fetal'])
            else:
                df.insert(0, 'geneset', [gene_set_names[gs.lower()] for gs in df['Annotation']])
        # Append to list
        genesets.append(df)

    # Concatenate gene lists if applicable
    if extract_geneset != 'all':
        geneset = genesets[0]
        # Return gene list or gene table
        if genes_only:
            return geneset['gene']
        else:
            return geneset
    else:
        # Reduce to gene name and geneset name
        for i in range(len(genesets)):
            genesets[i] = genesets[i][['geneset', 'gene']]
        # Return gene names only
        geneset = pd.concat(genesets)
        return geneset
    
    
def get_CanellasSocias_2022(species="Hs"):
    
    # Canellas-Socias, Nat, 2022
    # Data as listed in ST7 (only gene names)
    
    files = os.listdir(dir_gene_lists)
    if species == "Hs":
        reg = re.compile(rf'.*Canellas-Socias.*ST7\.csv') # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]), skiprows=1)
    else:
        if species == "Mm":
            reg = re.compile(rf'.*Canellas-Socias.*Mm\.csv') # Compile the regex
            lists = list(filter(reg.search, files))
            geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]))

            # Replace column names All HR to All-HR and TME HR to TME-HR
            geneset.columns = geneset.columns.str.replace('All.HR', 'All-HR')
            geneset.columns = geneset.columns.str.replace('TME.HR', 'TME-HR')
            # Replace all . by _
            geneset.columns = geneset.columns.str.replace('.', '_')
        else:
            raise AttributeError("Valid species are {Hs, Mm}.")
    
    return geneset

def get_collection_NV(species="Hs"):
    
    # Batlle-related gene sets as collected by Nuria V
    
    files = os.listdir(dir_gene_lists)
    if species == "Hs":
        reg = re.compile(rf'.*collection_NuriaV*.csv') # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]))
    else:
        raise AttributeError("Only valid species implemented: Hs")
    
    return geneset

def get_Centonze_2025(species="Hs"):

    # Centonze, Cancer Discov, 2025
    # Data as listed in ST1 (only gene names)
    # Gene sets in human or mouse nomenclature depending on primary source 
    
    files = os.listdir(dir_gene_lists)
    reg = re.compile(rf'.*Centonze.*ST1\.csv') # Compile the regex
    lists = list(filter(reg.search, files))
    if not lists:
        raise FileNotFoundError("Centonze gene set file not found.")
    geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]), skiprows=1, sep=';')
    
    return geneset

def get_Mustata_2013(directionality="up", species="Mm", genes_only=False):
    
    # Mustata, CellRep, 2013
    # Data as listed in ST1 with log2FC
    # (both up- and down-regulated genes available)
    
    # Valid input parameters
    if directionality not in ['up', 'down']:
        raise AttributeError("Valid directionalities are {up, down}.")
    if species not in ['Mm', 'Hs']:
        raise AttributeError("Valid species are {Mm, Hs}.")

    # Get list of all available files
    files = os.listdir(dir_gene_lists)
    
    # File identification
    if species == "Mm":
        reg = re.compile(rf'.*Mustata.*{directionality}\.csv') # Compile the regex
    elif species == "Hs":
         reg = re.compile(rf'.*Mustata.*{directionality}.*Hs\.csv') # Compile the regex
    lists = list(filter(reg.search, files))
    
    # Read mouse gene list and clean up
    if species == "Mm":
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]), skiprows=1)
        # Only keep columns with non-NA values
        geneset = geneset.iloc[:, :3]
        # Rename columns of dataframe
        colnames = {
            'Gene Symbol': 'gene',
            'Fold modulation (log2)  ': 'FC_log2',
            'Description': 'description',
        }
        geneset.rename(columns=colnames, inplace=True)
        
    # Read human gene list and clean up
    elif species == "Hs":
        geneset = pd.read_csv(os.path.join(dir_gene_lists, lists[0]))
        # Only keep columns with non-NA values        
        geneset = geneset.dropna(subset=['Human_symbol'])
    
    # Return gene list or gene table
    if genes_only:
        if species == "Mm":
            return geneset['gene']
        elif species == "Hs":
            return geneset['Human_symbol']
    else:
        return geneset


def get_Han_2018(data_id='DEGs', dir_data=None):
    
    # Han, Cell, 2018 [Mouse Cell Atlas]
    # Data as listed in ST4 or ST5
    
    # Set directory of gene lists if not provided
    if dir_data is None:
        dir_data = dir_gene_lists

    # Search for requested dataset
    files = os.listdir(dir_data)
    if (data_id == 'DEGs'):
        reg = re.compile(rf'.*Han.*DEGs\.csv')
    elif (data_id == 'SI'):
        reg = re.compile(rf'.*Han.*SI\.csv')
    else:
        raise ValueError("Data set not specified.")
    lists = list(filter(reg.search, files))    
    geneset = pd.read_csv(os.path.join(dir_data, lists[0]))

    # Index pd with gene names
    geneset.set_index('gene', drop=True, inplace=True)
    
    # Only keep columns with non-NA values
    geneset = remove_nan_cols(geneset)
    
    # Drop further obsolete columns
    if (data_id == 'DEGs'):
        geneset.drop(columns=['Unnamed: 0'], inplace=True)
    
    # Fill Cell Type column based on value above
    geneset['Cell Type'].fillna(method='ffill', inplace=True)
    
    return geneset


def get_Tirosh_2016(species="Hs"):
    
    # Tirosh, 2016, Nature
    # Data as listed in ST1 (only gene names)
    
    dir_gene_list = '/home/m014f/gene_lists/gene_lists_cell-cycle'
    files = os.listdir(dir_gene_list)
    
    if species == "Hs":
        reg = re.compile(rf'.*Tirosh.*ST1\.csv') # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_list, lists[0]))
    elif species == "Mm":
        reg = re.compile(rf'.*Tirosh.*ST1_convMH\.csv') # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_list, lists[0]))

        # Rename column G2_M to G2/M
        geneset.rename(columns={'G2_M': 'G2/M', 'G1_S': 'G1/S'}, inplace=True)
    else:
        raise AttributeError(f"Function not implemented yet for species {species}")
    # Only keep columns with non-NA values
    geneset = remove_nan_cols(geneset)
    
    # Clean-up leading spaces in each column
    #geneset = geneset.apply(strip_gene_names, axis=0)
    # Note: does not work yet, for now manually # TBC
    
    return geneset


def get_Hai_Hoffmann_2023(geneset, species="Hs"):
    
    # Hai, Hoffmann, Nat Commun, 2023
    # Data as listed in ST3-5
    
    dir_gene_list = '/home/m014f/gene_lists/gene_lists_brain/'
    files = os.listdir(dir_gene_list)
    
    genesets = {"Caprola6": "ST3",
                "CaprolaOn": "ST4",
                "Calcium": "ST5",}
    if geneset not in genesets.keys():
        print(f"Invalid geneset. Available options: {', '.join(genesets.keys())}")
    
    if species == "Hs":
        reg = re.compile(rf"HaiHoffmann.*{genesets[geneset]}.*\.csv") # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_list, lists[0]), skiprows=1)
    else:
        raise AttributeError(f"Function not implemented yet for species {species}")
    
    # Only keep columns with non-NA values
    geneset = remove_nan_cols(geneset)    
    
    return geneset


def get_RuizMoreno_2022(level = None, species="Hs"):
    
    # Hai, Hoffmann, Nat Commun, 2023
    # Data as listed in ST3-5
    
    dir_gene_list = '/home/m014f/gene_lists/gene_lists_brain/'
    files = os.listdir(dir_gene_list)
    
    if species == "Hs":
        reg = re.compile(rf"RuizMoreno.*\.csv") # Compile the regex
        lists = list(filter(reg.search, files))
        geneset = pd.read_csv(os.path.join(dir_gene_list, lists[0]), skiprows=1)
        
        # Rename columns
        cols = ['p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene']
        
        if level is None:
            levels = ['level_1', 'level_2', 'level_3', 'level_4']
            colnames = [f"{level}_{col}" for level in levels for col in cols]
            geneset.columns = colnames
        else:
            if level in [1,2,3,4]:
                geneset = geneset.iloc[:, ((level-1)*7):(level*7)]
                geneset.dropna(how='all', inplace=True)
                geneset.columns = cols
            else:
                print("Warning: Level value is not valid.")
                return
    else:
        raise AttributeError(f"Function not implemented yet for species {species}")
    
    # Only keep columns with non-NA values
    geneset = remove_nan_cols(geneset)    
    
    return geneset


def get_immunecell_scores(data="Fges_Charo"):
    
    dir_gene_lists = "/home/m014f/gene_lists/gene_lists_immune/"
    if data == "Fges_Charo":
        file = "Fges_Charo.gmt"
    else:
        if data == "Bindea":
            file = "Bindea.gmt"
        else:
            raise ValueError("Data set not valid.")
    if file not in os.listdir(dir_gene_lists):
        raise FileNotFoundError("Gene list not found.")

    # Read dataset
    print(f"Reading gene list from file {file}")
    gene_sets = []
    with open(os.path.join(dir_gene_lists, file), 'r') as f:
        for line in f:
            # Split each line into gene set name, description, and genes
            tokens = line.strip().split('\t')

            # Validate if there are enough tokens to avoid IndexError
            if len(tokens) < 3:
                print(f"Skipping line with insufficient data: {tokens}")
                continue  # Skip lines with insufficient elements
            
            # Extract gene set features
            gene_set = {
                'name': tokens[0],
                'class': tokens[1],
                'genes': tokens[2:]  # All remaining tokens are gene names
            }
            gene_sets.append(gene_set)
    
    # Return list of dictionaries
    print("Genesets provided as list of dictionaries.")
    return gene_sets


def get_TROP2KO_DGE(species="Mm", genotype="akps", genes_only = False, p_val=0.05):
    
    # TROP2 KO organoid bulk RNA-seq
    # obtained from Nuria

    dir_gene_lists = "/home/m014f/gene_lists/Jackstadt/E03/"

    # Valid input parameters
    if species not in ['Mm', 'Hs']:
        raise AttributeError("Valid species are {Mm, Hs}.")

    if genotype not in ['akps', 'kpn']:
        raise AttributeError("Valid genotypes are {akps, kpn}.")

    files = os.listdir(dir_gene_lists)

    # File identification based on species and genotype
    if species == "Mm":
        reg = re.compile(re.compile(rf'^tacstd2.*{genotype}.*\.tsv'))
        delim = "\t"
    elif species == "Hs":
        reg = re.compile(rf'^tacstd2.*{genotype}.*Hs\.csv')
        delim = ","
    files = list(filter(reg.search, files)) 
    
    # Reading the file
    geneset = pd.read_csv(os.path.join(dir_gene_lists, files[0]), delimiter=delim)

    # Filtering by adjusted p_value
    if p_val:
        print(f"The gene list is being filtered by adjusted p_value < {p_val}.")
        geneset[geneset['p_val_adj'] < p_val]
    else:
        print("Warning: gene list has not been filtered for statistically significant results.")

    # Sorting the gene set by log2FoldChange
    geneset = geneset.sort_values(by='avg_log2FC', ascending=False)

    '''
    if genes_only:
        if species == "Mm":
            genes = geneset['Mouse_symbol']
        else:
            genes = geneset['Human_symbol']
            # drop NA in gene list
            genes = genes[~genes.isna()]
            print(f"{len(genes)} of {len(geneset)} mouse genes have a known human homologue and are being returned.")
        return genes
    else:
        return geneset
    '''
    return geneset
    


def pre_network_checks(net, var_names, source, target, weight, label=None):
    ### Network checks
    if source not in net.columns and 'source' not in net.columns:
        raise ValueError(f"Source column {source} not found in network")
    if target not in net.columns and 'target' not in net.columns:
        raise ValueError(f"Target column {target} not found in network")
    if weight not in net.columns and 'weight' not in net.columns:
        warnings.warn(f"Weight column {weight} not found in network. Setting to None.")
        weight = None
    # Remove duplicate entries if available
    net = net[~net.duplicated([source, target])]
    # Check if gene list is non-zero
    if net.shape[0] == 0:
        warnings.warn(f"No pathways found in network {label}. Skipping...")
    # Check for each geneset if gene overlap is > 5. Otherwise exclude geneset.
    for geneset in net[source].unique():
        if len(set(net[net[source] == geneset][target]).intersection(var_names)) < 5:
            warnings.warn(f"Skipping {geneset} (less than 5 overlapping genes)")
            net = net[net[source] != geneset]
    # Return if no genesets are left
    if net.shape[0] == 0:
        raise ValueError(f"No pathways found in network {label}. Skipping...")
    else:
        # Rename columns to default if needed
        col_dict = {
            source: "source",
            target: "target",
            weight: "weight"
        }
        net = net.rename(columns=col_dict)
        return net
    

# Gene set names for publication
gs_names_for_pub = {
    # HR signatures from Canellas-Socias et al. Nat, 2022 @Batlle, ST7
    "score_allHR": "allHR (Cañellas-Socias et al.)",
    "score_coreHR": "coreHR (Cañellas-Socias et al.)",
    "score_epiHR": "epiHR (Cañellas-Socias et al.)",
    "score_tmeHR": "tmeHR (Cañellas-Socias et al.)",
    "score_Lgr5_Batlle": "Lgr5 (Cañellas-Socias et al.)",
    # Signatures listed in ST7 refrerring to original studies
    "score_Wnt_Batlle": "Wnt-On (Morral et al.)", # Morral et al. Cell Stem Cell, 2020 @Batlle
    "score_mKi67_Batlle": "mKi67 (Basak et al.)", # Basak et al. EMBO J, 2014 @Clevers
    "score_Basal-like_PDAC": "Basal-like pancreatic cancer (Raghavan et al.)", # Raghavan et al. Cell, 2021 @Shalek
    "score_Yap_Batlle": "YAP_22 (Wang et al.)", # Wang et al., Cell Rep, 2018 @Liang
       
    # Signatures from Álvarez-Varela et al. Nature Cancer, 2022 @Batlle, ST6 
    "score_Mucosecreeting": "Mucosecreeting (Álvarez-Varela et al.)",
    "score_Paneth_Cells": "Paneth cells (Álvarez-Varela et al.)",
    "score_Enteroendocrine": "Enteroendocrine (Álvarez-Varela et al.)",
    "score_Secretory_Progenitors": "Secretory progenitors (Álvarez-Varela et al.)",
    "score_Goblet_Cells": "Goblet cells (Álvarez-Varela et al.)",
    "score_Lgr5_signature": "Lgr5 (Álvarez-Varela et al.)",
    "score_LGR5Hi_MEX3AHi": "LGR5Hi MEX3AHi (Álvarez-Varela et al.)",
    "score_LGR5hi_MEX3Alow": "LGR5hi MEX3Alow (Álvarez-Varela et al.)",
    "YAP direct targets": "YAP targets (Álvarez-Varela et al.)",
    "score_YAP_direct_targets": "YAP direct targets (Álvarez-Varela et al.)",
    # Signatures listed in ST6 of Álvarez-Varela et al. Cell Stem Cell, 2022, referring to original studies
    # "score_Label_Retaining_Cells": "Label-retaining cells (Buczacki et al.)", # Buczacki et al., Nature, 2013 @Winton
    "score_Label_Retaining_Cells": "LRC (Buczacki et al.)", # Buczacki et al., Nature, 2013 @Winton
    "score_Crypt_proliferation": "Crypt proliferation (Jung et al.)", # Jung et al., Nature Medicine, 2011 @Batlle/Clevers 
    "score_Mex3A": "MEX3A (Barriga et al.)", # Barriga et al. Cell Stem Cell, 2017 @Batlle
    "score_mKI67_high_low": "mKi67 high-vs-low (Cortina et al.)", # Cortina et al. EMBO Mol Med, 2017 @Batlle

    # Plasticity signatures from Moorman et al. Nature, 2025 @Ganesh/Pe’er
    "score_fetal_Moorman": "Fetal (Moorman et al.)",
    "score_squamous_Moorman": "Squamous (Moorman et al.)",
    "score_neuroendocrine_Moorman": "Neuroendocrine (Moorman et al.)",

    # YAP KO and YAP OE signatures from Gregorieff et al. Nature, 2015 @Wrana
    "score_YAP": "YAP (Gregorieff et al.)",

    # Intestinal SC signatures from Munoz et al. EMBO Journal, 2012 @Clevers
    "score_Lgr5": "Lgr5 (Muñoz et al.)",

    # Fetal organoid signatures from Mustata et al. Cell Reports, 2013 @Garcia
    "score_fetal": "Fetal (Mustata et al.)",

    # Pol1R signature from Morral et al., Cell Stem Cell, 2019 @Batlle, ST5
    # Import from Nuria’s collection, better format.
    # get_collection_NV
    "score_POLR1A_High": "POLR1A High (Morral et al.)", # N=488

    
    # Implemented in a021.
    # EMP1 KRAS study from Centonze et al. Cancer Discovery, 2026 @Batlle, ST1
    # Signatures listed in ST1, referring to original study
    "score_iCMS2_(Joanito)": "iCMS2 (Joanito et al.)", # Joanito et al., Nature, 2022 @Tan
    "score_iCMS3_(Joanito)": "iCMS3 (Joanito et al.)", # Joanito et al., Nature, 2022 @Tan
    "score_Proliferation_(Merlos)": "Proliferation (Merlos-Suárez et al.)", # Merlos-Suárez et al., Cell Stem Cell, 2011 @Batlle
    "score_Immature_enterocytes_(Smillie)": "Immature enterocytes (Smillie et al.)", # Smillie et al., Cell, 2019 @Regev
    "score_revival_SCs_(Vazquez)": "Revival SCs (Vazquez et al.)", # Vazquez et al., Cell Stem Cell, 2022 @Leedham
    "score_Revival_SCs_(Vazquez)": "Revival SCs (Vazquez et al.)", # Vazquez et al., Cell Stem Cell, 2022 @Leedham
    "score_revival_SCs_(Ayyaz)": "Revival SCs SSC2c (Ayyaz et al.)", # Ayyaz et al., Nature, 2019 @Gregorieff
    "score_Revival_SCs_(Ayyaz)": "Revival SCs SSC2c (Ayyaz et al.)", # Ayyaz et al., Nature, 2019 @Gregorieff
    # Signatures from ST1, referring to original studies
    "score_EpiHR_(Cañellas)": 'EpiHR (Cañellas et al. 2022)',
    "score_Lgr5_ISCs_(Muñoz)": "Lgr5+ Intestinal Stem Cells (Muñoz et al. 2012)",
    "score_Fetal_Organoids_(Mustata)": "Fetal Organoids (Mustata et al. 2013)",
    "score_YAP_22_(Wang)": "YAP_22 (Wang et al. 2018)",
    "score_Basal_PDAC_(Raghavan)": "Basal pancreatic cancer cells (Raghavan et al. 2021)",
}