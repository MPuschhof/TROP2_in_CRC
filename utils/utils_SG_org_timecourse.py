import pandas as pd
# import genelist
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
import matplotlib.colors as mcolors

color_code = {
    'drug': {
        'order': ['UNADC', 'SG'],
        'colors': {
            'UNADC': '#9d9da1',
            'SG': '#ad7bb6',
        },
    },
}

def get_gene_sets(adata):

    import utils_genelist as genelist

    # Initiate gene set collection
    gene_sets = []

    ### YAP score (Gregorieff, Nature, 2015)
    print("Get Gregorieff YAP gene set")
    score_name = "score_YAP"
    gene_list = genelist.get_Gregorieff_2015(species="Hs", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

    ### Fetal score (Mustata, CellRep, 2013)
    print("Process fetal gene set")
    score_name = "score_fetal"
    gene_list = genelist.get_Mustata_2013(species="Hs", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

    ### Lgr5 score (Muñoz, EMBOJ, 2011)
    print("Process LGR5 gene set")
    score_name = "score_Lgr5"
    gene_list = genelist.get_Munoz_2011(species="Hs", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list)

    ### Ganesh-related gene sets
    print("Process Ganesh gene set collection")
    gene_list = genelist.get_Moorman_2023(extract_geneset="all", genes_only = True)
    for score_name in gene_list['geneset'].unique():
        df = gene_list[gene_list['geneset'] == score_name]
        genes = df['gene'].to_list()
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, genes)

    ### Emp1 score (Cañellas-Socias, Nature, 2022)
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
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list[col_name])

    ### Batlle-related gene set collection of Nuria V
    print("Process Batlle-related gene set collection of Nuria V")
    gene_list = genelist.get_collection_NV()
    # Iterate over all gene sets
    for score_name in gene_list.columns:
        gene_sets = genelist.add_to_gene_sets(adata, gene_sets, 'score_'+score_name, gene_list[score_name])    

    ### Flatten to one gene / row dataframe
    gene_sets = pd.DataFrame([
        {**d, 'gene': gene} for d in gene_sets for gene in d['genes']
    ]).drop('genes', axis=1)
    # Reorder columns
    gene_sets = gene_sets[['gene_set', 'gene', 'weight']]

    # Add results to anndata slot
    adata.uns['gene_sets'] = gene_sets

    return adata

# Custom gene set selection 
gene_set_selection = [
    # "score_allHR",
    "score_coreHR",
    "score_epiHR",
    # "score_tmeHR",
    "score_Lgr5_Batlle",
    "score_Wnt_Batlle",
    "score_mKi67_Batlle",
    # "score_Basal-like_PDAC", # listed again below.
    "score_Yap_Batlle", # listed again below.
    "score_Mucosecreeting",
    "score_Paneth_Cells",
    "score_Enteroendocrine",
    "score_Secretory_Progenitors",
    "score_Goblet_Cells",
    "score_Lgr5_signature", 
    "score_LGR5Hi_MEX3AHi",
    "score_LGR5hi_MEX3Alow",
    "score_YAP_direct_targets",
    "score_fetal_Moorman",
    "score_squamous_Moorman",
    "score_neuroendocrine_Moorman",
    "score_YAP",
    "score_Lgr5",
    "score_fetal",
    "score_POLR1A_High",
    # "score_iCMS2_(Joanito)",
    # "score_iCMS3_(Joanito)",
    "score_Proliferation_(Merlos)",
    "score_Immature_enterocytes_(Smillie)",
    "score_Revival_SCs_(Vazquez)",
    "score_Revival_SCs_(Ayyaz)",
    # "score_EpiHR_(Cañellas)", already listed above.
    # "score_Lgr5_ISCs_(Muñoz)", already listed above.
    # "score_Fetal_Organoids_(Mustata)",  already listed above.
    "score_YAP_22_(Wang)",
    "score_Basal_PDAC_(Raghavan)",    
]


# Gene set names for publication
gs_names_for_pub = {
    # HR signatures from Canellas-Socias et al. Nat, 2022 @Batlle, ST7
    "score_allHR": "allHR (Cañellas-Socias et al.)",
    "score_coreHR": "coreHR (Cañellas-Socias et al.)",
    "score_epiHR": "epiHR (Cañellas-Socias et al.)",
    "score_tmeHR": "tmeHR (Cañellas-Socias et al.)",
    "score_Lgr5_Batlle": "Lgr5 (Cañellas-Socias et al.)", ### only 5 genes
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
    "score_Goblet_Cells": "Goblet cells (Álvarez-Varela et al.)", ### only 5 genes
    "score_Lgr5_signature": "Lgr5 (Álvarez-Varela et al.)", ### only 5 genes
    "score_LGR5Hi_MEX3AHi": "LGR5Hi MEX3AHi (Álvarez-Varela et al.)",
    "score_LGR5hi_MEX3Alow": "LGR5hi MEX3Alow (Álvarez-Varela et al.)",
    "YAP direct targets": "YAP targets (Álvarez-Varela et al.)", ### only 7 genes
    "score_YAP_direct_targets": "YAP direct targets (Álvarez-Varela et al.)",
    # Signatures listed in ST6 of Álvarez-Varela et al. Cell Stem Cell, 2022, referring to original studies
    "score_Label_Retaining_Cells": "Label-retaining cells (Buczacki et al.)", # Buczacki et al., Nature, 2013 @Winton
    "score_Crypt_proliferation": "Crypt proliferation (Jung et al.)", # Jung et al., Nature Medicine, 2011 @Batlle/Clevers 
    "score_Mex3a_Barriga": "MEX3A (Barriga et al.)", # Barriga et al. Cell Stem Cell, 2017 @Batlle
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
    "score_Revival_SCs_(Vazquez)": "Revival SCs (Vazquez et al.)", # Vazquez et al., Cell Stem Cell, 2022 @Leedham
    "score_Revival_SCs_(Ayyaz)": "Revival SCs SSC2c (Ayyaz et al.)", # Ayyaz et al., Nature, 2019 @Gregorieff 
    # Signatures from ST1, referring to original studies
    "score_EpiHR_(Cañellas)": 'EpiHR (Cañellas et al. 2022)',
    "score_Lgr5_ISCs_(Muñoz)": "Lgr5+ Intestinal Stem Cells (Muñoz et al. 2012)",
    "score_Fetal_Organoids_(Mustata)": "Fetal Organoids (Mustata et al. 2013)",
    "score_YAP_22_(Wang)": "YAP_22 (Wang et al. 2018)",
    "score_Basal_PDAC_(Raghavan)": "Basal pancreatic cancer cells (Raghavan et al. 2021)",
}