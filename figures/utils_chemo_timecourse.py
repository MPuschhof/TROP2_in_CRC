import os
import pandas as pd
import scanpy as sc
from pathlib import Path
# import genelist
import matplotlib.pyplot as plt
import numpy as np

import utils
import utils_genelist as genelist

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
    'treatment': {
        'order': ['drug', 'drug_cont', 'drug_release'],
        'colors': {
            'drug': ylorrd[3],
            'drug_cont': "#8e44ad",
            'drug_release': blues[1],
        } 
    },
    'treatment_binary': { # 2026-01-12
        'order': ['treatment', 'release'],
        'colors': {
            # get related but distinct colors: darker than sandybrown
            'treatment': 'peru',
            'release': 'peachpuff',
        } 
    },
    'treatment_categ': {
        'order': ['none', 'chemo', 'release'],
        'colors': {
            'none': 'lightgrey',
            'chemo': 'sandybrown',
            'release': 'peachpuff',
        }   
    },
    'time': {
        'order': order_time,
        'colors': dict(zip(order_time, colors_HTO_maxID))
    },
}


def add_time(adata):
    adata.obs['time'] = [hto.rsplit("-",1)[1] for hto in adata.obs['HTO_maxID']]
    adata.obs['time'] = pd.Categorical(adata.obs['time'],
                                            categories = order_time,
                                            ordered = True)
