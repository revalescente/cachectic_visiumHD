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
from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf


sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26foxO.zarr")

sdata['blocco4_intissue'] = sdata['blocco4_intissue'][sdata['blocco4_intissue'].name=="c26foxO"]

sopa.aggregate(sdata, key_added = 'nuclei_counts_nop', bins_key= "filtered", 
shapes_key = "blocco4_nuclei_boundaries", expand_radius_ratio=1, min_transcripts=1, 
min_intensity_ratio=0.1, no_overlap = True)
