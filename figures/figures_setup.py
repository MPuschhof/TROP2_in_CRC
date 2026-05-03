# Packages used throughout
import scanpy as sc
import anndata
import numpy as np
import scipy as sp
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import os
import sys
from pathlib import Path
import session_info

# Setting data path
dir_data = Path(f"{os.getcwd()}/../data_arrayExpress")


# # Local modules used throughout
# module_path = os.path.abspath('../src')
# if module_path not in sys.path:
#     sys.path.append(module_path)