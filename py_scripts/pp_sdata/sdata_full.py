import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
from spatialdata import SpatialData
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import sopa
import re
import json
from shapely.affinity import scale
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized
import py_scripts.pp_sdata.pp_functions as pp

# read data with io_reader
# block_numbers = [1, 2, 3, 4, 5, 6, 7, 9]
block_numbers = [4, 7, 9]
spe_blocks = {}

for i in block_numbers:
    block_name = f"spe_blocco{i}"
    spe_blocks[block_name] = visium_hd(
        path=f"/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/blocco{i}/outs",
        dataset_id=f"blocco{i}",
        filtered_counts_file=False,
        bin_size='002',
        bins_as_squares=True,
        annotate_table_by_labels=False,
        fullres_image_file=f"/mnt/europa/valerio/HE_images/color_corrected/pp_blocco{i}_20x.tif",
        load_all_images=False,
        var_names_make_unique=True
    )

# Write data as zarr stores, spe_blocks is a dictionary of datasets
for block_name, spe in spe_blocks.items():
    spe.write(f"/mnt/europa/valerio/data/{block_name}.zarr")

# apply the sdata preprocess function to filter in tissue bins and save it back to a zarr file
for num in block_numbers:
    path = f'/mnt/europa/valerio/data/zarr_store/general/spe_blocco{num}.zarr'
    try:
        result = pp.sdata_pp(path)
        # (Optional) Save or use result here
        print(f"Processed blocco {num}")
        result.write(f'/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco{num}.zarr')
        print(f"Zarr Store created for blocco{num}")
    except Exception as e:
        print(f"Error processing blocco {num}: {e}")

# read it back in a list 
block_numbers = [1, 2, 3, 4, 5, 6, 7, 9]
#block_numbers = [4, 7, 9]
spe_blocks_list = []
for i in block_numbers:
  spe_blocks_list.append(sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco{i}.zarr'))

# or a dictionary
# block_numbers = [4, 7, 9]
spe_blocks = {}
for i in block_numbers:
  block_name = f"blocco{i}"
  spe_blocks[block_name] = sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco{i}.zarr')

# ------------------------------------------------------------------------------

# Blocks merging ---------------------------------------------------------------
# Don't know if it will be useful
# spatialdata.concatenate(sdatas, region_key=None, instance_key=None, concatenate_tables=False, 
#                         obs_names_make_unique=True, modify_tables_inplace=False, attrs_merge=None)
# join the blocks
# spe = sd.concatenate(spe_blocks_list, concatenate_tables = True) # se l'oggetto è una lista non rinomina nulla
# spe.write('/mnt/europa/valerio/data/cachetic_fulldataset.zarr')

# Single Sample objects --------------------------------------------------------
# separate the blocks' samples
# SHAM:                  sham
# C26:                   c26
# C26+shMuRF1/shpanFoxO: c26murf1
# C26+shpanFoxO:         c26foxO
# C26+shSMAD2/3:         c26SMAD23
# C26+shSTAT3:           c26STAT3


# Read the dictionary with the sample-blocks-bbox information

with open('/mnt/europa/valerio/data/blocco_sample_bbox_dict.json', 'r') as f:
    blocco_sample_bbox_dict = json.load(f)

subset_dict = {name: blocco_sample_bbox_dict[name] for name in ['blocco4', 'blocco7', 'blocco9'] if name in blocco_sample_bbox_dict}
# for all dataset use: blocco_sample_bbox_dict
pp.divide_samples(spe_blocks, subset_dict)

# problem: blocco7 don't find sham sample after the bbox filtering - SOLVED (sham and c26 inverted in dict)
# spe_blocks['blocco7']['blocco7_intissue_002um']['name'] = spe_blocks['blocco7']['blocco7_intissue_002um']['name'].astype('category')
# plt.figure(figsize=(50, 50))
# ax = plt.gca()
# spe_blocks['blocco7'].pl.render_shapes('blocco7_intissue_002um', color = "name", outline=False, outline_alpha=0, outline_width=0, fill_alpha=1).pl.show(ax = ax, save = "output_python/roba_da_buttare/is_sham_inb7.png")


# we have the zarr store for each sample of the blocks. Now check if it's working

spe_blocks = {}
# change with blocco_sample_bbox_dict for full dataset 
for blocco, samples in subset_dict.items():
    for sample in samples:
        zarr_name = f"{blocco}_{sample}"
        zarr_path = f"/mnt/europa/valerio/data/zarr_store/blocchi/{zarr_name}.zarr"
        spe_blocks[zarr_name] = read_zarr_standardized(zarr_path)

read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/.zarr")

# ok, creating and writing this data work, check with a plot

# rasterize --------------------------------------------------------------------
# rasterize bins to fast rapresentation
# rasterize_bins() requires a compresed sparse column (csc) matrix
# spatialdata.rasterize_bins(sdata, bins, table_name, col_key, row_key, value_key=None, return_region_as_labels=False)
# for blocco, samples in subset_dict.items():
#   for sample in samples:
#     sample_name = f"{blocco}_{sample}"
#     spe_blocks[sample_name]["filtered"].X = spe_blocks[sample_name]["filtered"].X.tocsc()
#     spe_blocks[sample_name]["rasterized_bins"] = sd.rasterize_bins(
#         spe_blocks[sample_name],
#         f"{blocco}_intissue_002um",
#         "filtered",
#         "array_col",
#         "array_row",
#     )
# 
# sample_names = [
#     "blocco4_c26SMAD23", "blocco4_c26foxO", "blocco4_c26",
#     "blocco7_sham", "blocco7_c26STAT3", "blocco7_c26",
#     "blocco9_c26murf1", "blocco9_c26foxO", "blocco9_c26SMAD23"
# ]
# 
# for sample_name in sample_names:
#     blocco = sample_name.split('_')[0]
#     fig, ax = plt.subplots(figsize=(15, 15))
#     spe_blocks[sample_name].pl.render_images(
#         "rasterized_bins", channel="Trim63", scale="full"
#     ).pl.show(coordinate_systems=blocco, ax=ax, 
#     save = f"output_python/roba_da_buttare/{sample_name}_trim63.png"
#     )
#     plt.close(fig)  # Close figure to avoid memory leaks
# 
# # problem: i didn't remove the bins belonging to other exp_cond
# plt.figure(figsize=(50, 50))
# ax = plt.gca()
# spe_blocks['blocco4_c26foxO'].pl.render_shapes("blocco4_intissue_002um", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
# ).pl.show(ax = ax, coordinate_systems = "blocco4", save="output_python/roba_da_buttare/b4_c26foxO.png")
# 
# 
# # let's keep only the bins inside the right tissue
# spe_blocks['blocco4_c26foxO']['blocco4_intissue_002um'] = spe_blocks['blocco4_c26foxO']['blocco4_intissue_002um'][spe_blocks['blocco4_c26foxO']['blocco4_intissue_002um']['name'] == 'c26foxO']
# spe_blocks['blocco4_c26foxO']['filtered'] = sd.match_table_to_element(spe_blocks['blocco4_c26foxO'], element_name="blocco4_intissue_002um", table_name="filtered")







