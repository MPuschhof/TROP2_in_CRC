import os
import pandas as pd
import numpy as np
from pathlib import Path
import utils
import utils_genelist as genelist

### Preprocessing parameters
# QC metric: mitochondrial reads
mt_thr = 30 # set to 30 for all after obtaining S03 results
# Number of genes / cell
gene_qnt = 0.95 # upper threshold
min_genes = 750 # lower threshold
# Scrublet doublet score
doublet_thr = .15
# Genes picked up in at least min_cells cells
min_cells = 3 # lower threshold


### C01-a006
def get_meta_data_C01a006(adata, sample_name):
    adata.obs['sample_id'] = sample_name
    adata.obs['donor_id'] = sample_name.split("_")[0]
    # Define and extract experiment information
    exp_dict = {
        "SF1.1": ["HD4246_1336179", "HD4246_1337158", "HD4246_1336139",],
        "SF3.0": ["HD4254_1336213", "HD4254_1337131", "HD4254_1337186",],
        "NV28": ["HD4277_1277095","HD4277_1277808",],
        "NV41": ["HD4309_1327464", "HD4309_1327524", "HD4309_1327521", "HD4309_1327518", "HD4309_1327519", "HD4309_1327520",],
        "SF4.0": ["HD4362_1355736"],
        "NV21": ["HD4362_1272219", "HD4362_1272215",],
        "NV40": ["HD4236_1327502", "HD4236_1327505", "HD4236_1327494", "HD4236_1327504", "HD4236_1327506",],
        "MM4.1.1": ["HD4272_1350128", "HD4272_1350215", "HD4272_1350216", "HD4272_1350217",],
    }    
    adata.obs['experiment_id'] = [utils.find_key_for_value(exp_dict, sample) for sample in adata.obs['sample_id']]

    return adata
# def get_meta_data_C01a006(bc_data, sample_name):
#     bc_data['sample_id'] = sample_name
#     bc_data['donor_id'] = sample_name.split("_")[0]
#     # Define and extract experiment information
#     exp_dict = {
#         "SF1.1": ["HD4246_1336179", "HD4246_1337158", "HD4246_1336139",],
#         "SF3.0": ["HD4254_1336213", "HD4254_1337131", "HD4254_1337186",],
#         "NV28": ["HD4277_1277095","HD4277_1277808",],
#         "NV41": ["HD4309_1327464", "HD4309_1327524", "HD4309_1327521", "HD4309_1327518", "HD4309_1327519", "HD4309_1327520",],
#         "SF4.0": ["HD4362_1355736"],
#         "NV21": ["HD4362_1272219", "HD4362_1272215",],
#         "NV40": ["HD4236_1327502", "HD4236_1327505", "HD4236_1327494", "HD4236_1327504", "HD4236_1327506",],
#         "MM4.1.1": ["HD4272_1350128", "HD4272_1350215", "HD4272_1350216", "HD4272_1350217",],
#     }    
#     bc_data['experiment_id'] = [utils.find_key_for_value(exp_dict, sample) for sample in bc_data['sample_id']]

#     return bc_data

# # annos_full_dataset = {
# #     'filename': 'Biobank_subcut_14_merged_harmony',
# #     'resolution': 0.3,
# #     'labels': {
# #         '0': 'Tumor cells 1',
# #         '1': 'Immune cells',
# #         '2': 'Tumor cells 2',
# #         '3': 'Tumor cells 3',
# #         '4': 'Tumor cells 4',
# #         '5': 'Tumor cells 5',
# #         '6': 'Fibroblasts',
# #         '7': 'Tumor cells 6',
# #         '8': 'Endothelial cells',
# #     }
# # }

# annos_full_dataset = {
#     'filename': 'Biobank_subcut_24_merged_harmony',
#     'resolution': 0.3,
#     'labels': {
#         '0': 'Tumor cells 1',
#         '1': 'Tumor cells 2',
#         '2': 'Tumor cells 3',
#         '3': 'Immune cells',
#         '4': 'Tumor cells 4',
#         '5': 'Tumor cells 5',
#         '6': 'Fibroblasts',
#         '7': 'Tumor cells 6',
#     }
# }

# annos_epi_dataset = {
#     'filename': 'Biobank_subcut_24_epi_harmony',
#     'resolution': 0.4,
#     'labels': {
#         '0': 'LGR5+', # Lgr5, Wnt-Batlle
#         '1': 'revSC', # CLU # Mucosecreting, LGR5Hi-Mex3AHi
#         '2': 'LGR5+ ribosomal', # none.
#         '3': 'TROP2+, KRT20+', #HRCs, YAP fetal, LGR5Hi-Mex3AHi
#         '4': 'TROP2+, KRT20-', # HRCs, YAP
#         '5': 'Prolif. LGR5+', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A, LRCs, Lgr5
#         '6': 'Prolif.', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A
#         '7': 'MEX3A+', # Lgr5Hi-Mex3AHi,
#     }
# }

# # annos_epi_dataset = { # 241119-21
# #     'filename': 'Biobank_subcut_24_epi_harmony',
# #     'resolution': 0.4,
# #     'labels': {
# #         '0': 'LGR5+ (1)', # Lgr5, Wnt-Batlle
# #         '1': '(Secr prog.)', # Mucosecreting, LGR5Hi-Mex3AHi
# #         '2': 'LGR5+ (2)', # none.
# #         '3': 'TROP2+ (1), KRT20+', #HRCs, YAP fetal, LGR5Hi-Mex3AHi
# #         '4': 'TROP2+ (2), KRT20-', # HRCs, YAP
# #         '5': 'Prolif. (1), LGR5+', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A, LRCs, Lgr5
# #         '6': 'Prolif. (2)', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A
# #         '7': 'MEX3A+', # Lgr5Hi-Mex3AHi,
# #     }
# # }

# # annos_epi_dataset = { # v01
# #     'filename': 'Biobank_subcut_24_epi_harmony',
# #     'resolution': 0.4,
# #     'labels': {
# #         '0': 'LGR5+ (1)', # Lgr5, Wnt-Batlle
# #         '1': '(Secr prog.)', # Mucosecreting, LGR5Hi-Mex3AHi
# #         '2': 'LGR5+ (2)', # none.
# #         '3': 'TROP2+ HRCs (1), KRT20+', #HRCs, YAP fetal, LGR5Hi-Mex3AHi
# #         '4': 'TROP2+ HRCs (2), KRT20-', # HRCs, YAP
# #         '5': 'Prolif. (1), LGR5+', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A, LRCs, Lgr5
# #         '6': 'Prolif. (2)', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A
# #         '7': 'NE / MEX3A+', # Lgr5Hi-Mex3AHi,
# #     }
# # }

# annos_epi_dataset_Batlle_style = { # v02
#     'filename': 'Biobank_subcut_24_epi_harmony',
#     'resolution': 0.4,
#     'labels': {
#         '0': 'LGR5+ KI67- (1)', # Lgr5, Wnt-Batlle
#         '1': '(Secr prog.)', # Mucosecreting, LGR5Hi-Mex3AHi
#         '2': 'LGR5+ KI67- (3)', # none.
#         '3': 'HRCs KRT20+', # HRCs, YAP fetal, LGR5Hi-Mex3AHi
#         '4': 'HRCs KRT20-', # HRCs, YAP
#         '5': 'KI67+ (1)', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A, LRCs, Lgr5
#         '6': 'KI67+ (2)', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A
#         '7': 'Undef.', # Lgr5Hi-Mex3AHi,
#     }
#     # Comparison LGR5 clusters (0 vs 2):
#     # Cl 0 enriched: WNT5B, ZNRF3, LGR5, WNT2B, NOTUM, LRIG1, TIAM1, TCF7L2, ...
#     # Cl 2 enriched: REG4, LYZ, AGR2, LY6E, TFF1, ribosomal genes, ...

#     # Comparison TROP2 clusters (3 vs 4):
#     # Cl 3 enriched: ITGAM, EGR3, FOSB, ATF3, TEAD2,  ...
#     # Cl 4 enriched: ITGA5, FABP1, LRP1/4, CA9, ALDOA, AGR2, TFF1/3, LRIG2, CEACAM1, LYZ ...
# }

# # annos_epi_dataset_res03 = {
# #     'filename': 'Biobank_subcut_14_merged_harmony_subset_epi_harmony',
# #     'resolution': 0.3,
# #     'labels': {
# #         '0': 'LGR5+', # heavy in RPS/RPLs
# #         # indications of Notch (DLL4, Hes1) and ATF4 (JUN)
# #         '1': 'TROP2+ (KRT20+)',
# #         '2': 'LGR5+ (KRT20+)',
# #         '3': 'KI-67+',
# #         '4': 'Diff.', # score LGR5+ MEX3A+
# #     }
# # }

# # annos_epi_dataset_res05 = {
# #     'filename': 'Biobank_subcut_14_merged_harmony_subset_epi_harmony',
# #     'resolution': 0.5,
# #     'labels': {
# #         '0': 'LGR5+ (1)',
# #         '1': 'LGR5+ (2)',
# #         '2': 'TROP2+ KRT20+ HRCs',
# #         '3': 'TROP2+',
# #         '4': 'TROP2+ KI-67+', # score LGR5+ MEX3A-low
# #         '5': 'KI-67+',
# #         '6': 'Diff.', # score LGR5+ MEX3A+
# #     }
# # }

# def get_celltypist_model(species="Hs"):
#     # # Get celltypist model for immune cells
#     if species == 'Mm':
#         ct_models = {
#             "m_ad_gut": "/home/m014f/references/celltypist/Adult_Mouse_Gut.pkl",
#         }
#     elif species == 'Hs':
#         ct_models = {
#             "h_cells_int": "/home/m014f/references/celltypist/Cells_Intestinal_Tract.pkl",
#             "h_ad_int": "/home/m014f/references/celltypist/Adult_Human_Intestine.pkl",
#         }
#     else:
#         raise ValueError("Species not implemented.")

#     # Check if files exist
#     for model_path in ct_models.values():
#         if not Path(model_path).is_file():
#             raise FileNotFoundError(f"Model file not found: {model_path}")

#     print(f"Intestinal celltypist models (species '{species}') loaded.")
#     return ct_models

# def get_gene_sets(adata):

#     # Initiate gene set collection
#     gene_sets = []

#     ### YAP score (Gregorieff, Nature, 2015)
#     print("Get Gregorieff YAP gene set")
#     score_name = "score_YAP"
#     gene_list = genelist.get_Gregorieff_2015(species="Hs", genes_only = True)
#     gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

#     ### Fetal score (Mustata, CellRep, 2013)
#     print("Process fetal gene set")
#     score_name = "score_fetal"
#     gene_list = genelist.get_Mustata_2013(species="Hs", genes_only = True)
#     gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

#     ### Lgr5 score (Muñoz, EMBOJ, 2011)
#     print("Process LGR5 gene set")
#     score_name = "score_Lgr5"
#     gene_list = genelist.get_Munoz_2011(species="Hs", genes_only = True)
#     gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

#     ### Ganesh-related gene sets
#     print("Process Ganesh gene set collection")
#     gene_list = genelist.get_Moorman_2023(extract_geneset="all", genes_only = True)
#     for score_name in gene_list['geneset'].unique():
#         df = gene_list[gene_list['geneset'] == score_name]
#         genes = df['gene'].to_list()
#         gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes)

#     ### Emp1 score (Cañellas-Socias, Nature, 2022)
#     print("Process EMP1-related gene sets")
#     scores = {
#         "score_allHR": 'All-HR',
#         "score_coreHR": 'coreHRC',
#         "score_epiHR": 'EpiHR',
#         "score_tmeHR": 'TME-HR',
#         "score_Lgr5_Batlle": 'Lgr5 signature',	
#         "score_Wnt_Batlle": 'Wnt ON signature (Morral et al. 2020)',
#         "score_mKi67_Batlle": 'mKi67 (Basak et al. 2014)',
#         "score_Yap_Batlle": 'YAP_22 (Wang et al. 2018)'
#     }
#     gene_list = genelist.get_CanellasSocias_2022(species="Hs")
#     # Iterate over all scores
#     for score_name,col_name in scores.items():
#         gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list[col_name])

#     ### TROP2KO bulk RNA_seq score (Jackstadt lab)
#     # Get TROP2KO genes for both genotypes
#     genes_akps = genelist.get_TROP2KO_DGE(species="Hs", genotype = "akps", p_val=0.05)
#     genes_kpn = genelist.get_TROP2KO_DGE(species="Hs", genotype = "kpn", p_val=0.05)
#     # Derive up- and downregulated genes upon TROP2 KO
#     genes_akps_KO = genes_akps[(genes_akps['avg_log2FC'] < 0)]
#     genes_akps_wt = genes_akps[(genes_akps['avg_log2FC'] > 0)]
#     genes_kpn_KO = genes_kpn[(genes_kpn['avg_log2FC'] < 0)]
#     genes_kpn_wt = genes_kpn[(genes_kpn['avg_log2FC'] > 0)]
#     # Derive union of marker genes
#     genes_KO = list(set(genes_akps_KO['Human_symbol']).union(set(genes_kpn_KO['Human_symbol'])))
#     genes_wt = list(set(genes_akps_wt['Human_symbol']).union(set(genes_kpn_wt['Human_symbol'])))
#     # Calculate gene set scores
#     score_name = "score_TROP2KO_up"   
#     gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes_KO)
#     score_name = "score_TROP2KO_dn"   
#     gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes_wt)

#     ### Batlle-related gene set collection of Nuria V
#     print("Process Batlle-related gene set collection of Nuria V")
#     gene_list = genelist.get_collection_NV()
#     # Iterate over all gene sets
#     for score_name in gene_list.columns:
#         gene_sets = genelist.add_to_gene_sets(adata, gene_sets, 'score_'+score_name, gene_list[score_name])    

#     ### Flatten to one gene / row dataframe
#     gene_sets = pd.DataFrame([
#         {**d, 'gene': gene} for d in gene_sets for gene in d['genes']
#     ]).drop('genes', axis=1)
#     # Reorder columns
#     gene_sets = gene_sets[['gene_set', 'gene', 'weight']]

#     # Add results to anndata slot
#     adata.uns['gene_sets'] = gene_sets

#     return adata


# def add_SC_status(adata):
#     # Define gene sets and colors
#     features = ["LGR5", "TACSTD2"]
#     colors = ['green', "red"]
#     # For each cell check if feature gene expression is > 0
#     for feature, color in zip(features, colors):
#         if feature in adata.var_names:
#             adata.obs[f"{feature}_pos"] = adata.obs_vector(feature, layer = "counts") > 0
#         else:
#             print(f"Feature {feature} not found in dataset.")
#     ### Assign SC population status using gene expression levels
#     conditions = [
#         (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == True),
#         (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == False),
#         (adata.obs[f"{features[0]}_pos"] == False) & (adata.obs[f"{features[1]}_pos"] == True),
#     ]
#     # Corresponding choices for each condition
#     choices = ["both", features[0], features[1]]
#     # Create the new column based on the conditions and choices
#     adata.obs['SC_status'] = np.select(conditions, choices, default="none")
#     return adata


# def add_TROP2_status(adata):
#     # Define gene sets and colors
#     feature = "TACSTD2"
#     adata.obs['TROP2_status'] = ['TROP2_pos' if val > 0 else 'TROP2_neg' for val in adata.obs_vector(feature, layer = "counts")]
#     return adata

# def add_LGR5_status(adata):
#     # Define gene sets and colors
#     feature = "LGR5"
#     adata.obs['LGR5_status'] = ['LGR5_pos' if val > 0 else 'LGR5_neg' for val in adata.obs_vector(feature, layer = "counts")]
#     return adata

# def add_EMP1_status(adata):
#     # Define gene sets and colors
#     feature = "EMP1"
#     adata.obs['EMP1_status'] = ['EMP1_pos' if val > 0 else 'EMP1_neg' for val in adata.obs_vector(feature, layer = "counts")]
#     return adata

# def add_TROP2_EMP1_status(adata):
#     # Define gene sets and colors
#     features = ["EMP1", "TACSTD2"]
#     colors = ['blue', "green"]
#     # For each cell check if feature gene expression is > 0
#     for feature, color in zip(features, colors):
#         if feature in adata.var_names:
#             adata.obs[f"{feature}_pos"] = adata.obs_vector(feature, layer = "counts") > 0
#         else:
#             print(f"Feature {feature} not found in dataset.")
#     ### Assign SC population status using gene expression levels
#     conditions = [
#         (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == True),
#         (adata.obs[f"{features[0]}_pos"] == True) & (adata.obs[f"{features[1]}_pos"] == False),
#         (adata.obs[f"{features[0]}_pos"] == False) & (adata.obs[f"{features[1]}_pos"] == True),
#     ]
#     # Corresponding choices for each condition
#     choices = ["both", features[0], features[1]]
#     # Create the new column based on the conditions and choices
#     adata.obs['TROP2_EMP1_status'] = np.select(conditions, choices, default="none")
#     return adata