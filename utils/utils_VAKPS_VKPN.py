import pandas as pd
import utils_genelist as genelist

anno_leiden_dict_epi = {
    # 'filename': 'seuratObj_MP_epi_res0.6',
    'resolution': 0.6,
    'labels': {
        '0': 'ISC', # Lgr5+
        '1': 'Emp1+ Krt20+', # or HRC KRT20+?
        '2': 'Mitochondrial',
        '3': 'ISC/HRC', # Lgr5+ Tacstd2+
        '4': 'Fetal/HRC KRT20-', # Tacstd2+ Emp1+
        '5': 'Prolif. HRC', # Prolif.
        '6': 'Prolif. ISC', # Prolif.
        '7': 'Wnt Ascl2+', 
    }
}

color_code = {
    'SC_status': {
        'order': ['Lgr5', 'both',  'none', 'Tacstd2'],
        'colors': ['red', 'orange', 'grey', 'green'],
    },
    'TROP2_status': {
        'order': ['positive', 'negative'],
        'colors': {
            'positive': 'green',
            'negative': 'red'
        }
    }
}

# Get gene sets
def get_gene_sets(adata, verbose=True):

    # Initiate gene set collection
    gene_sets = []

    ### Fetal score (Mustata, CellRep, 2013)
    if verbose:
        print("Process fetal gene set")
    score_name = "score_fetal"
    gene_list = genelist.get_Mustata_2013(species="Mm", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list, verbose=verbose)

    ### Lgr5 score (Muñoz, EMBOJ, 2011)
    if verbose:
        print("Process LGR5 gene set")
    score_name = "score_Lgr5"
    gene_list = genelist.get_Munoz_2011(species="Mm", genes_only = True)
    gene_sets = genelist.add_to_gene_sets(adata, gene_sets, score_name, gene_list, verbose=verbose)


    ### Flatten to one gene / row dataframe
    gene_sets = pd.DataFrame([
        {**d, 'gene': gene} for d in gene_sets for gene in d['genes']
    ]).drop('genes', axis=1)
    # Reorder columns
    gene_sets = gene_sets[['gene_set', 'gene', 'weight']]

    # Add results to anndata slot
    adata.uns['gene_sets'] = gene_sets

    return adata
