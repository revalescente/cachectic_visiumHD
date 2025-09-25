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
import sys

# --- MODIFICATION ---
# The full path is now passed directly from the bash script,
# so we can use it directly instead of constructing it.
if len(sys.argv) < 2:
    print("Error: Please provide the full path to the .zarr file.")
    sys.exit(1)

path_sdata = sys.argv[1]
# path_sdata = "/mnt/europa/valerio/data/zarr_store/blocchi/*.zarr"

print(f"Processing file: {path_sdata}")

sdata = read_zarr_standardized(path_sdata)
sdata = sf.postprocess_step(sdata)

print(f"Finished processing: {path_sdata}")
