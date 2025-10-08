import spatialdata as sd
from spatialdata_io import visium_hd
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
import geopandas as gpd
import pandas as pd
import os
import sopa
import re
import json
from shapely.affinity import scale
from sopa.io.standardize import sanity_check, read_zarr_standardized
import py_scripts.pp_sdata.pp_functions as pp
from py_scripts.utils.utils_fun import read_from_json

# 1. read data and filtering, then saving it back in a new zarr
# select the number of blocks to pp if you want
pp.preprocess_step(all = False)

# 2. dividing samples into individual zarr stores

# let's use a dictionary to keep track of every sample for every blocco 
dict_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json"
samples_dict = read_from_json(dict_path)

# if you wanna use it for a subset
subset_dict = {name: blocco_sample_bbox_dict[name] for name in ['blocco7', 'blocco9'] if name in blocco_sample_bbox_dict}
# for all dataset use: blocco_sample_bbox_dict
pp.divide_samples(samples_dict)


# ------------------------------------------------------------------------------

# if you want to use all the samples at once

# spe_blocks = {}
# # change with blocco_sample_bbox_dict for full dataset 
# for blocco, samples in subset_dict.items():
#     for sample in samples:
#         zarr_name = f"{blocco}_{sample}"
#         zarr_path = f"/mnt/europa/valerio/data/zarr_store/blocchi/{zarr_name}.zarr"
#         spe_blocks[zarr_name] = read_zarr_standardized(zarr_path)
# 
# 

