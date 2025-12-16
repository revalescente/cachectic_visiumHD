import spatialdata as sd
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import (Identity, set_transformation)
import numpy as np
import geopandas as gpd
import pandas as pd
import os
#import gc
import sopa
import json
#import time
from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
# from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
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
# path_sdata = "/mnt/europa/valerio/data/zarr_store/samples/*.zarr"

print(f"Processing file: {path_sdata}")

sdata = sd.read_zarr(path_sdata)
#sdata = sf.segmentation_step(sdata)
sdata = sf.postprocess_step(sdata, expand_radius_ratio = 0, min_transcripts = 4, min_intensity_ratio = 0.15, no_overlap = True)

filtered_nuclei_key = "filtered_nuclei"

sdata["nuclei_counts"].obs["region"] = filtered_nuclei_key
sdata.set_table_annotates_spatialelement("nuclei_counts", region = filtered_nuclei_key, region_key="region", instance_key="cell_id")
  

sdata['nuclei_counts'] = sd.match_table_to_element(sdata, element_name = filtered_nuclei_key, table_name='nuclei_counts')
  

sdata['intissue_002um']['cell_id'] = sf.bin_to_cell_id_vector(sdata, table_key = 'nuclei_counts')
sdata['filtered_bins'] = sdata['intissue_002um'][sdata['intissue_002um']['cell_id'].notna()]

sdata.write_element('filtered_bins')
sdata.write_element('filtered_nuclei', overwrite = True)
sdata.write_element('nuclei_counts', overwrite = True)
sdata.write_element('raster_nuclei', overwrite = True)

sdata.delete_element_from_disk('image_patches')

print(f"Finished processing: {path_sdata}")
