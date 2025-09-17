import spatialdata as sd
import spatialdata_plot
from spatialdata import SpatialData
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import (Identity, set_transformation)
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import sopa
import json
from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf
new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
# from importlib import reload
# reload(sf)

# read dict
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    blocco_sample_bbox_dict = json.load(f)

for blocco, samples_dict in blocco_sample_bbox_dict.items():
    # Now iterate through each sample in this blocco
    for sample_key in samples_dict.keys():
        # Construct the path for this specific sample
        path_sdata = f"/mnt/europa/valerio/data/zarr_store/blocchi/{blocco}_{sample_key}.zarr"
        sdata = read_zarr_standardized(path_sdata)
        sdata = sf.segmentation_step(sdata)
  
# -------------------------------------------------------------------------------
# Testing - to cancel all 
# sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26SMAD23.zarr")



for blocco, samples_dict in blocco_sample_bbox_dict.items():
    # Now iterate through each sample in this blocco
    for sample_key in samples_dict.keys():
        # Construct the path for this specific sample
        path_sdata = f"/mnt/europa/valerio/data/zarr_store/blocchi/{blocco}_{sample_key}.zarr"
        print(path_sdata)





