import anndata
import spatialdata as sd
import numpy as np
import geopandas as gpd
import pandas as pd
import anndata as ad
import os
import sopa
import re
import sys
import json
import shutil
from collections import Counter
# Automatically add repo root to path to fix "No module named 'py_scripts'"
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if repo_root not in sys.path:
    sys.path.append(repo_root)

from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp
import py_scripts.segmentation.segm_functions as sf

# 1. Get the blocco key from the command line
if len(sys.argv) < 2:
    print("Usage: python process_single_blocco.py <BLOCCO_KEY>")
    sys.exit(1)

BLOCCO_KEY = sys.argv[1]

# disable autosave of sopa
sopa.settings.auto_save_on_disk = False

# paths of interest
sdata_dir = "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins"
output_dir1 = "/mnt/europa/valerio/data/zarr_store/adatas/binned_adatas/version_2.0.0/016um"
output_dir2 = "/mnt/europa/valerio/data/zarr_store/adatas/binned_adatas/version_2.0.0/008um"


# Load dictionary
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

if BLOCCO_KEY not in samples_dict:
    print(f"[ERROR] {BLOCCO_KEY} not found in dictionary.")
    sys.exit(1)

block_samples = samples_dict[BLOCCO_KEY]

print(f"\n{'='*50}")
print(f"=== PROCESSING BLOCK: {BLOCCO_KEY} ===")
print(f"{'='*50}")

# table_name = "arivis_nuclei_table"
table_name1 = "square_016um"
table_name2 = "square_008um"

for sample_name, sample_info in block_samples.items():
    print(f"\n--- Processing sample: {BLOCCO_KEY} - {sample_name} ---")
    sdata = sd.read_zarr(f"{sdata_dir}/{BLOCCO_KEY}_{sample_name}")
    # add centroid coords in the .obs
    # sdata[table_name].obs['y_coord'] = sdata[table_name].obsm['spatial'][:, 0]
    # sdata[table_name].obs['x_coord'] = sdata[table_name].obsm['spatial'][:, 1]
    # remove and resave the table with the last mods
    # sdata.delete_element_from_disk(table_name)
    # sdata.write_element(table_name)
    #export as adata to ease
    adata1 = sdata[table_name1].copy()
    adata1.obs['y_coord'] = adata1.obsm['spatial'][:, 0]
    adata1.obs['x_coord'] = adata1.obsm['spatial'][:, 1]
    adata2 = sdata[table_name2].copy()
    adata2.obs['y_coord'] = adata2.obsm['spatial'][:, 0]
    adata2.obs['x_coord'] = adata2.obsm['spatial'][:, 1]
    del adata1.obsm
    del adata2.obsm
    adata1.write_h5ad(f'{output_dir1}/{BLOCCO_KEY}_{sample_name}.h5ad')
    adata2.write_h5ad(f'{output_dir2}/{BLOCCO_KEY}_{sample_name}.h5ad')

# Must remove obsm otherwise anndataR as_sce doesn't work... sadly
