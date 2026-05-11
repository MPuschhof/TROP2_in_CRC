anno_leiden_dict = {
    # 'filename': 'HD4246_merged_harmony_res04',
    'resolution': 0.4,
    'labels': {
        '0': 'revSC', #'CLU+',
        '1': 'Fetal', #TACSTD2+ (fetal-like)',
        '2': 'Prolif. ISC', #'LGR5+ prolif.',
        '3':  'HRC', # 'EMP1+ (non-prolif.)',
        '4': 'Ribosomal',
        '5': 'Prolif.',
        '6': 'ISC', # 'LGR5+ (non-prolif.)',
        '7': 'Fetal/HRC', # 'TACSTD2+ EMP1+ (non-prolif.)',
        '8': 'EMT', # 'Undef.',
    }
}

color_code = {
    'treatment_resp': {
        'order': ['down', 'up', 'stable',],
        'colors': ['#5B84B1FF', '#FC766AFF', 'grey']
    },
    'SC_status': {
        'order': ['LGR5', 'both',  'none', 'TACSTD2'],
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


# Custom gene set selection 
custom_subset = [
    "score_coreHR",
    "score_epiHR",
    "score_Lgr5_Batlle",
    "score_Wnt_Batlle",
    "score_mKi67_Batlle",
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
    "score_Label_Retaining_Cells",
    "score_Proliferation_(Merlos)",
    "score_Immature_enterocytes_(Smillie)",
    "score_Revival_SCs_(Vazquez)",
    "score_Revival_SCs_(Ayyaz)",
    "score_YAP_22_(Wang)",
    "score_Basal_PDAC_(Raghavan)",    
]