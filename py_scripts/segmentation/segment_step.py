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
import gc
import sopa
import json
import time
from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf
new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
# from importlib import reload
# reload(sf)

# parellelization settings
sopa.settings.parallelization_backend = "dask"

# For processing large H&E images
sopa.settings.dask_client_kwargs = {
    "n_workers": 2,                    # Fewer workers for larger memory per worker
    "memory_limit": "20GB",            # More memory per worker
    "worker_options": {
        "memory_target_fraction": 0.7,  # Spill to disk at 70% memory usage
        "memory_spill_fraction": 0.8,   # Pause worker at 80% usage
        "memory_pause_fraction": 0.95   # Terminate tasks at 95% usage
    }
}

# read dict
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    blocco_sample_bbox_dict = json.load(f)

# subset_dict = {name: blocco_sample_bbox_dict[name] for name in ['blocco4', 'blocco5'] if name in blocco_sample_bbox_dict}

for blocco, samples_dict in subset_dict.items():
    # Now iterate through each sample in this blocco
    for sample_key in samples_dict.keys():
        # Construct the path for this specific sample
        path_sdata = f"/mnt/europa/valerio/data/zarr_store/blocchi/{blocco}_{sample_key}.zarr"
        sdata = read_zarr_standardized(path_sdata)
        sdata = sf.segmentation_step(sdata)
        
        # Explicitly delete and collect garbage
        del sdata
        gc.collect()
        time.sleep(5)
        
# individual sample
path_sdata = "/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_sham.zarr"
sdata = read_zarr_standardized(path_sdata)
sdata = sf.segmentation_step(sdata)
# sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26SMAD23.zarr")






