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


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

dir_figures = Path(__file__).resolve().parent
dir_repo = dir_figures.parent

dir_data = dir_repo / "data_arrayExpress"
dir_utils = dir_repo / "utils"

# Make repo-level utils folder importable
if str(dir_utils) not in sys.path:
    sys.path.append(str(dir_utils))