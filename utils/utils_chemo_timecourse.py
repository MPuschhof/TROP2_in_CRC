import pandas as pd
# import genelist
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
import matplotlib.colors as mcolors


# Experimental time points representing treatment options
shared = ['0h', '12h', '24h', '48h', '72h']
treatment = ['120h']
release = ['72h+12h', '72h+24h', '72h+48h', '72h+72h']
order_time = [tp for tp in shared + treatment + release]

### Generate color code for HTO_maxID
cmap = plt.get_cmap('YlOrRd')
ylorrd = [mcolors.to_hex(cmap(x)) for x in np.linspace(0.2, 0.9, 5)]
purple = "#8e44ad" 
cmap = plt.get_cmap('Blues')
blues = [mcolors.to_hex(cmap(x)) for x in np.linspace(0.9, 0.4, 4)]

# Combine into final color list
colors_HTO_maxID = ylorrd + [purple] + blues

color_code = {
    'time': {
        'order': order_time,
        'colors': dict(zip(order_time, colors_HTO_maxID))
    },
}


def add_time(adata):
    adata.obs['time'] = [hto.rsplit("-",1)[1] for hto in adata.obs['HTO_classification']]
    adata.obs['time'] = pd.Categorical(adata.obs['time'],
                                            categories = order_time,
                                            ordered = True)




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
    "score_allHR",
    "score_coreHR",
    "score_epiHR",
    "score_tmeHR",
    "score_Lgr5_Batlle",
    "score_Wnt_Batlle",
    "score_mKi67_Batlle",
    "score_Yap_Batlle",
    "score_Mucosecreeting",
    "score_Paneth_Cells",
    "score_Enteroendocrine",
    "score_Lgr5_signature", 
    "score_LGR5Hi_MEX3AHi",
    "score_LGR5hi_MEX3Alow",
    "score_fetal_Moorman",
    "score_squamous_Moorman",
    "score_neuroendocrine_Moorman",
    "score_YAP",
    "score_Lgr5",
    "score_fetal",
    "score_POLR1A_High", 
    "score_Label_Retaining_Cells",
    "score_Crypt_proliferation",
]

# Gene set names for publication
gs_names_for_pub = {
    # HR signatures from Canellas-Socias et al. Nat, 2022 @Batlle, ST7
    # get_CanellasSocias_2022
    "score_allHR": "All-HRC",
    "score_coreHR": "Core-HRC",
    "score_epiHR": "Epi-HRC",
    "score_tmeHR": "TME-HRC",
    "score_Lgr5_Batlle": "Lgr5 score (Batlle)",
    # Signatures listed in ST7 refrerring to original studies
    "score_Wnt_Batlle": "Wnt score (Batlle)", # Morral et al. Cell Stem Cell, 2020 @Batlle
    "score_mKi67_Batlle": "Ki67 (Battle)", # Basak et al. EMBO J, 2014 @Clevers
    "score_Yap_Batlle": "YAP (Batlle)", # Wang et al., Cell Rep, 2018 @Liang
       
    # # Signatures from Álvarez-Varela et al. Nature Cancer, 2022 @Batlle, ST6 
    "score_LGR5Hi_MEX3AHi": "LGR5+ MEX3A+",
    "score_LGR5hi_MEX3Alow": "LGR5+ MEX3A-low",
    "score_Crypt_proliferation": "Crypt proliferation", # Jung et al., Nature Medicine, 2011 @Batlle/Clevers 

    # Plasticity signatures from Moorman et al. Nature, 2025 @Ganesh/Pe’er
    # get_Moorman_2023
    "score_neuroendocrine_Moorman": "Neuroendocrine",

    # YAP KO and YAP OE signatures from Gregorieff et al. Nature, 2015 @Wrana
    # get_Gregorieff_2015
    "score_YAP": "YAP", 

    # Intestinal SC signatures from Munoz et al. EMBO Journal, 2012 @Clevers
    # get_Munoz_2011
    "score_Lgr5": "Lgr5 score", 

    # Fetal organoid signatures from Mustata et al. Cell Reports, 2013 @Garcia
    # get_Mustata_2013
    "score_fetal": "Fetal", 

    # Pol1R signature from Morral et al., Cell Stem Cell, 2019 @Batlle, ST5
    # get_collection_NV
    "score_POLR1A_High": "POLR1A+", # N=488
}