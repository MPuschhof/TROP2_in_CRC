anno_leiden_dict_full = { # Full data set
    # 'filename': 'Biobank_subcut_24_merged_harmony',
    'resolution': 0.3,
    'labels': {
        '0': 'Tumor cells 1',
        '1': 'Tumor cells 2',
        '2': 'Tumor cells 3',
        '3': 'Immune cells',
        '4': 'Tumor cells 4',
        '5': 'Tumor cells 5',
        '6': 'Fibroblasts',
        '7': 'Tumor cells 6',
    }
}

anno_leiden_dict_epi = { # Epithelial subset
    # 'filename': 'Biobank_subcut_24_epi_harmony',
    'resolution': 0.4,
    'labels': {
        '0': 'ISC', # Lgr5, Wnt-Batlle
        '1': 'revSC', # CLU # Mucosecreting, LGR5Hi-Mex3AHi
        '2': 'Ribosomal ISC', # Lgr5 ribosomal, 
        '3': 'Fetal/HRC', # 'TROP2+, KRT20+', #HRCs, YAP fetal, LGR5Hi-Mex3AHi
        '4': 'HRC', # 'TROP2+, KRT20-', # HRCs, YAP
        '5': 'Prolif. ISC', # 'Prolif. LGR5+', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A, LRCs, Lgr5
        '6': 'Prolif.', # 'Prolif.', # Lgr5Hi-Mex3ALow, Ki67, prolif, POLR1A
        '7': 'DTP',# 'MEX3A+', # Lgr5Hi-Mex3AHi,
    }
}

color_code = {
    'SC_status': {
        'order': ['LGR5', 'both',  'none', 'TACSTD2'],
        'colors': ['red', 'orange', 'grey', 'green'],
    },
    'EMP1_status': {
        'order': ['TACSTD2', 'both',  'none', 'EMP1'],
        'colors': ['green', 'orange', 'grey', 'royalblue'],
    }
}

def add_TROP2_EMP1_status(adata):
    # Define gene sets and colors
    features = ["EMP1", "TACSTD2"]
    colors = ['blue', "green"]
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
    adata.obs['TROP2_EMP1_status'] = np.select(conditions, choices, default="none")
    return adata


# Custom gene set selection 
custom_subset = [
    "score_YAP",
    "score_coreHR",
    "score_epiHR",
    "score_fetal",
    "score_LGR5Hi_MEX3AHi",

    "score_LGR5hi_MEX3Alow",
    "score_mKi67_Batlle",
    "score_squamous_Moorman",
    "score_Yap_Batlle",
    "score_POLR1A_High",

    "score_fetal_Moorman",
    "score_Paneth_Cells",
    "score_Crypt_proliferation",
    "score_Mucosecreeting",
    "score_neuroendocrine_Moorman",

    "score_Lgr5_Batlle",
    "score_Label_Retaining_Cells",
    "score_Wnt_Batlle",
    "score_Lgr5",
]