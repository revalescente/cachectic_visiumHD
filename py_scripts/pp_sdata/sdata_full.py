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
import shutil
from shapely.affinity import scale
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from py_scripts.pp_sdata.pp_functions import sdata_pp


# let's display the areas where no expression is detected as transparent
new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
new_cmap


# read data with visium reader
# spe_blocco1 = visium_hd(path = '/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/blocco1/outs', 
#                 dataset_id='blocco1', 
#                 filtered_counts_file=False, 
#                 bin_size='002', 
#                 bins_as_squares=True, 
#                 annotate_table_by_labels=False, 
#                 fullres_image_file='/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco1_20x.tif', 
#                 load_all_images=False, 
#                 var_names_make_unique=True)

block_numbers = [1, 2, 3, 4, 5, 6, 7, 9]
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
        fullres_image_file=f"/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco{i}_20x.tif",
        load_all_images=False,
        var_names_make_unique=True
    )

# Assuming spe_blocks is your dictionary of datasets
for block_name, spe in spe_blocks.items():
    spe.write(f"/mnt/europa/valerio/data/{block_name}.zarr")

for num in block_numbers:
    path = f'/mnt/europa/valerio/data/spe_blocco{num}.zarr'
    try:
        result = sdata_pp(path)
        # (Optional) Save or use result here
        print(f"Processed blocco {num}")
        # Delete the old zarr store
        result.write(f'/mnt/europa/valerio/data/filtered_blocco{num}.zarr')
        print(f"Zarr Store created for blocco{num}")
        # Delete the old zarr store (handles file or folder)
        if os.path.isdir(path):
            shutil.rmtree(path)
            print(f"Deleted directory {path}")
        elif os.path.isfile(path):
            os.remove(path)
            print(f"Deleted file {path}")
        else:
            print(f"Warning: {path} does not exist as file or directory")
    except Exception as e:
        print(f"Error processing blocco {num}: {e}")
        # Do NOT delete the file if there was an error

# read zarr filtered data
for block_name, spe in spe_blocks.items():
    spe.write(f"/mnt/europa/valerio/data/{block_name}.zarr")

spe_blocks = {}
for i in block_numbers:
  block_name = f"spe_blocco{i}"
  spe_blocks[block_name] = sd.read_zarr(f'/mnt/europa/valerio/data/filtered_blocco{i}.zarr')


















