import anndata
import spatialdata as sd
import numpy as np
import geopandas as gpd
import pandas as pd
import anndata as ad
import os
import sopa
import re
import json
import shutil
from collections import Counter
from py_scripts.utils.utils_fun import read_from_json


samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
samples_key = [
    details['sample_key'] 
    for blocco in samples_dict.values() 
    for details in blocco.values()
]
input_dir = "/mnt/europa/valerio/data/zarr_store/binned_samples/version_1.0.0"
output_dir = "/mnt/europa/valerio/data/zarr_store/adatas/binned_adatas"

# sdata = sd.read_zarr(f"{input_dir}/blocco1_sham")
# 
# sdata['filtered'].obs['y_coord'] = sdata['filtered'].obsm['spatial'][:, 0]
# sdata['filtered'].obs['x_coord'] = sdata['filtered'].obsm['spatial'][:, 1]

for sample in samples_key:
    # read sdata
    sdata = sd.read_zarr(f"{input_dir}/{sample}")
    # add centroid coords in the .obs
    sdata['filtered'].obs['y_coord'] = sdata['filtered'].obsm['spatial'][:, 0]
    sdata['filtered'].obs['x_coord'] = sdata['filtered'].obsm['spatial'][:, 1]
    # remove and resave the table with the last mods
    sdata.delete_element_from_disk('filtered')
    sdata.write_element('filtered')
    #export as adata to ease
    adata = sdata['filtered'].copy()
    adata.write_h5ad(f'{output_dir}/{sample}')
