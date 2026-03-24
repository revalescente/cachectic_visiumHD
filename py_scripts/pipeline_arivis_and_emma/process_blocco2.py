import sys
import os
import gc
import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import re
import numpy as np
import geopandas as gpd
import pandas as pd
from skimage import io
from shapely.affinity import scale
from spatialdata.models import ShapesModel, Image2DModel, Labels2DModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from skimage.measure import regionprops_table
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp
import py_scripts.segmentation.segm_functions as sf

BLOCCO_KEY = 'blocco2'

# disable autosave of sopa
sopa.settings.auto_save_on_disk = False

# paths of interest
intissue_gfp_dir = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
data_b2 = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell"
arivis_dir = "/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks"
fullres_dir = "/mnt/europa/valerio/HE_images/color_corrected/blocchi"
fluo_path = "/mnt/europa/valerio/Fluo_images/warped_tif/samples"
save_dir = "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins"

# Load dictionary
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

if BLOCCO_KEY not in samples_dict:
    print(f"[ERROR] {BLOCCO_KEY} not found in dictionary.")

block_samples = samples_dict[BLOCCO_KEY]

fullres_path = f"{fullres_dir}/pp_{BLOCCO_KEY}_20x.tif"
visium_path = f"{data_b2}/{BLOCCO_KEY}/outs"

sdata = visium_hd(
    path=visium_path,
    fullres_image_file=fullres_path,
    dataset_id=BLOCCO_KEY,
    filtered_counts_file=False,
    bin_size=['002','008','016'],
    bins_as_squares=True,
    annotate_table_by_labels=False,
    load_all_images=False,
    var_names_make_unique=True,
    image_models_kwargs={'dims' : ['c', 'y', 'x']},
    load_segmentations_only=False,
    load_nucleus_segmentations=False
)

# Remove the cell segmentations of spaceranger to save memory
if f'{BLOCCO_KEY}_cell_segmentations' in sdata:
    del sdata[f'{BLOCCO_KEY}_cell_segmentations']
if 'cell_segmentations' in sdata:
    del sdata['cell_segmentations']

# plottino

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('blocco2_full_image', scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = 'output_python/blocco2_new_pain.png')

# plot of the 3 tissue to check if they are corrects
# new version
{'c26': {'min_coordinate': [0, 13000], 'max_coordinate': [21414, 20069], 'sample_key': 'blocco2_c26'}, 
'c26murf1': {'min_coordinate': [0, 8750], 'max_coordinate': [21414, 13750], 'sample_key': 'blocco2_c26murf1'}, 
'c26SMAD23': {'min_coordinate': [0, 0], 'max_coordinate': [21414, 8750], 'sample_key': 'blocco2_c26SMAD23'}}

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 13000],
    max_coordinate=[21414, 20069],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p1.png")

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 8750],
    max_coordinate=[21414, 13750],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p2.png")

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 0],
    max_coordinate=[21414, 8750],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p3.png")

# old version
     "blocco2": {
         "c26": {
-            "min_coordinate": [0, 0],
-            "max_coordinate": [21414, 7200],
             "sample_key": "blocco2_c26"
         },
         "c26murf1": {
-            "min_coordinate": [0, 6500],
-            "max_coordinate": [21414, 12000],
             "sample_key": "blocco2_c26murf1"
         },
         "c26SMAD23": {
-            "min_coordinate": [0, 11500],
-            "max_coordinate": [21414, 20069],
             "sample_key": "blocco2_c26SMAD23"
         }
     },

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 0],
    max_coordinate=[21414, 7200],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p1_old.png")

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 6500],
    max_coordinate=[21414, 12000],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p2_old.png")

plt.figure(figsize=(10, 10))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 11500],
    max_coordinate=[21414, 20069],
    target_coordinate_system=BLOCCO_KEY
).pl.render_images("blocco2_full_image", scale = "scale2").pl.show(ax = ax, coordinate_systems=BLOCCO_KEY, save = "output_python/blocco2_p3_old.png")

# file json aggiornato per far tornare i tagli sull'immagine raddrizzata, ma io poi ero tornato a lavorare sull'immagine storta.. quindi dovevo solo recuperare le coordinate per
# i tagli della versione vecchia del file json... GRAZIE GIT!! 
# restore the old version of json 
