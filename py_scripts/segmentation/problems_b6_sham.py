import spatialdata as sd
#import spatialdata_plot
#from spatialdata import SpatialData
#from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
#from spatialdata.transformations import (Identity, set_transformation)
import numpy as np
#import scanpy as sc
#import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
#import gc
import sopa
import json
#import time
#from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
#from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
#from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf

sopa.settings.auto_save_on_disk = False

sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr")

sf.sopa_attrs_check(sdata)
sopa.make_image_patches(sdata, patch_width = 300, roi_key = "blocco6_intissue")
sopa.segmentation.stardist(sdata, model_type='2D_versatile_he', 
  image_key= sdata.attrs['cell_segmentation_image'], min_area=10, delete_cache=True, 
  recover=False, prob_thresh=0.2, nms_thresh=0.6, key_added = f'blocco6_nuclei_boundaries')

# let's plot all the important stuff

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images("blocco6_full_image", scale = "scale2").pl.render_shapes("blocco6_intissue", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.render_shapes("image_patches", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco6", save = 'output_python/b6_probs.png')


# let's do all from the filtered zarr

spe = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco6.zarr")
# "sham": {
#             "min_coordinate": [0, 5500],
#             "max_coordinate": [16144, 12000],
#             "sample_key": "blocco6_sham"
            
sdata_bbox = spe.query.bounding_box(
                axes=["x", "y"],
                min_coordinate= [0, 5500],
                max_coordinate= [16144, 12000],
                target_coordinate_system="blocco6"
)
# Subset elements
subset_keys = ["blocco6_full_image", "blocco6_intissue_002um", "blocco6_intissue", "filtered"]
sdata_subset = sdata_bbox.subset(subset_keys, filter_tables=False)

# filtering element not of my sample
sdata_subset["blocco6_intissue_002um"] = sdata_subset["blocco6_intissue_002um"][sdata_subset["blocco6_intissue_002um"]['name'] == "sham"]
sdata_subset["blocco6_intissue"] = sdata_subset["blocco6_intissue"][sdata_subset["blocco6_intissue"]['name'] == "sham"]

# matching table with shapes
sdata_subset["filtered"] = sd.match_table_to_element(sdata_subset, element_name="blocco6_intissue_002um", table_name="filtered")

sopa.utils.set_sopa_attrs(
                sdata_subset,
                cell_segmentation_key = "blocco6_full_image",
                tissue_segmentation_key = "blocco6_full_image",
                bins_table_key = "filtered"
)

sanity_check(sdata_subset)

sdata_subset.write("/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr")

