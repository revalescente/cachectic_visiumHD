import anndata
import spatialdata as sd
import numpy as np
import geopandas as gpd
import pandas as pd
import os
import sopa
import re
import json
import shutil
from collections import Counter
from py_scripts.utils.utils_fun import read_from_json

# dizionario campioni
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

# extract the sample_key in a list
samples_key = [
    details['sample_key'] 
    for blocco in samples_dict.values() 
    for details in blocco.values()
]

# leggiamo come lista
# samples = ["blocco1_c26foxO", "blocco1_sham", "blocco3_c26murf1"]
spe_blocks_list = []

for sample in samples:
    spe_blocks_list.append(sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/samples/{sample}.zarr'))

# leggiamo come dizionario 
# subset the sdata with the only table we are interested in, then concat

# elements names
el_list = ['full_image', 'raster_nuclei', 'filtered_bins', 'filtered_nuclei', 'nuclei_counts']

spe_blocks = {}
for sample in samples:
  sdata = sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/samples/{sample}.zarr')
  subdata = sdata.subset(el_list)
  spe_blocks[sample] = subdata

# spe_dict = sd.concatenate(spe_blocks, concatenate_tables = True, region_key = "region", instance_key = "cell_id")

# concateniamo i due e vediamo le differenze
spe_list = sd.concatenate(spe_blocks_list, concatenate_tables = True)
spe_dict = sd.concatenate(spe_blocks, 
    region_key = "region", 
    instance_key = "cell_id",
    concatenate_tables = True 
)

spe_list.write('/mnt/europa/valerio/data/zarr_store/concat_all_samples_almost.zarr')

# some problem to concatenate object with a dictionary (all object have the same name/objects, should be straight forward)

# ------------------------------------------------------------------------
# Extract the 'nuclei_counts' table from each sample
nuclei_counts_list = [
    spe_blocks[sample]['nuclei_counts']
    for sample in samples
]

# Single samples
for adata, sample in zip(nuclei_counts_list, samples):
  adata.obsm.clear()
  pathname = f"/mnt/europa/valerio/data/zarr_store/adatas/{sample}.h5ad"
  adata.write(pathname)

for i, adata in enumerate(anndata_list):
    adata.obsm.clear()  # This removes all elements from .obsm
    fname = f"sample_{i+1}.h5ad"
    adata.write(fname)

# Concatenate all AnnData objects into one
nuclei_counts_concat = anndata.concat(nuclei_counts_list, 
                                      join="outer",  # or "inner" if you want only shared genes
                                      label="sample_id", 
                                      keys=samples)

# nuclei_counts_concat is now a single AnnData object with all nuclei_counts
print(nuclei_counts_concat)
nuclei_counts_concat.write("/mnt/europa/valerio/data/zarr_store/concat_allsamples.h5ad")

# cannot read in R if there are spatial elements, so remove "spatial" 
# from obsm and add in the sce as a spatial experiment / spatial feature experiment

spatial_coords = nuclei_counts_concat.obsm.pop('spatial')

for key in ['spatial', 'bins_assignments', 'intensities']:
    if key in nuclei_counts_concat.obsm:
        del nuclei_counts_concat.obsm[key]

# as a single cell experiment
nuclei_counts_concat.write("/mnt/europa/valerio/data/zarr_store/concat_as_sce.h5ad")

# with spatial coordinate on a separate object
np.savetxt("spatial_coords.csv", spatial_coords, delimiter=",")










