# Leiden annotations full dataset
anno_leiden_full = {
    # 'filename': 'Biobank09_crg-adj_merged',
    'resolution': 0.6,
    'labels': {
        '0': 'T cells',
        '1': 'Monocytes / Macrophages',
        '2': 'Epithelial cells',
        '3': 'T cells',
        '4': 'T cells',
        '5': 'B cells',
        '6': 'Plasma cells',
        '7': 'T cells',
        '8': 'T cells',
        '9': 'Pericytes',
        '10': 'T cells',
        '11': 'Fibroblasts',
        '12': 'Plasmacytoid dendritic cells',
        '13': 'Endothelial cells',
        '14': 'Mast cells',
        '15': 'Monocytes / Macrophages',
        '16': 'Enteric glia cells',
    }
}

markers_leiden_full = {
    "T cells": ["CD3D", "TRAC", "CD8B"], 
    "Monocytes / Macrophages": ["CD14", "S100A9", "LYZ"], 
    "Epithelial cells": ["EPCAM", "KRT8", "KRT19"],
    "B cells": ["MS4A1", "CD19", "VPREB3"], 
    "Plasma cells": ["SDC1", "CD38", "CD27"], 
    "Pericytes": ["PDGFRB", "RGS5", "CSPG4"], 
    "Fibroblasts": ["PDGFRA", "FAP", "THY1"],
    "Plasmacytoid dendritic cells": ['IRF7', 'GZMB', 'GPR183'], 
    "Endothelial cells": ["PECAM1", "VWF", "KDR"],
    "Mast cells": ["KIT", "TPSAB1", "CPA3"],
    "Enteric glia cells": ["GFRA3", "S100B"],
}


color_code = {
    'SC_status': {
        'order': ['LGR5', 'both',  'none', 'TACSTD2'],
        'colors': ['red', 'orange', 'grey', 'green'],
    }
}

# Order of donors for Fig S1N
order_donor_id = [
    'HD4254',
    'HD4181',
    'HD4249',
    'HD4309',
    'HD4241',
    'HD4362',
    'HD4246',
    'HD4272',
]
