import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
#import geopandas as gpd
#import pandas as pd
import sopa
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized

# Set parallelization backend to Dask
# sopa.settings.parallelization_backend = 'dask'     # not working

# read filtered data
filtered = sd.read_zarr(os.path.expanduser("~/data/b1_cleaned.zarr")) # non mi chiedere perché 

# setting up sopa metadata
sopa.utils.set_sopa_attrs(
    filtered,
    cell_segmentation_key="blocco1_full_image",
    tissue_segmentation_key="blocco1_full_image",
    bins_table_key="cleaned",
    boundaries_key = "blocco1_intissue"
)
sanity_check(filtered)
# metadata for single channel stardist segmentation
# sopa.utils.set_sopa_attrs(
#     filtered,
#     cell_segmentation_key="blocco1_conv_image",
#     tissue_segmentation_key="blocco1_full_image",
#     bins_table_key="cleaned"
# )


# subset the data so we keep only one tissue for this test
# y = (0,8000), x as is

# subsetting for testing
b1_stat3 = filtered.query.bounding_box(
    min_coordinate=[0, 0],
    max_coordinate=[16166, 8000],
    axes=("x", "y"),
    target_coordinate_system="blocco1",
)


# Check if a SpatialData object meets Sopa's requirements
sanity_check(b1_stat3)

# Write a standardized SpatialData object to disk
write_standardized(b1_stat3, os.path.expanduser("~/data/b1_stat3.zarr"))

#-------------------------------------------------------------------------------

# read standardized object with sopa function
b1_stat3 = read_zarr_standardized(os.path.expanduser("~/data/b1_stat3.zarr"))

# segmentation 

# divide the image in patches to segments individually first
# sopa.make_image_patches(sdata, patch_width=2000, patch_overlap=50, image_key=None, key_added=None)
sopa.make_image_patches(b1_stat3, patch_width = 500, patch_overlap = 50)

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images('blocco1_conv_image', cmap = "grey"
).pl.render_shapes("image_patches", fill_alpha = 0, outline_alpha = 1).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_convolved.png')


# let's see how sopa tissue segmentation work
sopa.segmentation.tissue(b1_stat3, expand_radius_ratio=0.05, key_added = "sopa_intissue")

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "sopa_intissue", outline_alpha=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_intissue.png')

# my tissue seg
b1_stat3['sub_intissue'] = b1_stat3['blocco1_intissue'].iloc[[0]]

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "sub_intissue", outline_alpha=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_myintissue.png')

# not working, but we already have the intissue shapes

# stardist segmentation of the convolved image
sopa.segmentation.stardist(b1_stat3, model_type='2D_versatile_fluo', channels = "conv", image_key= 'blocco1_conv_image', 
min_area=10, delete_cache=True, recover=False, prob_thresh=0.2, nms_thresh=0.6)

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images('blocco1_full_image').pl.render_shapes('stardist_boundaries', fill_alpha = 0, outline_alpha = 1
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_segments_bad.png')

# not even close

# stardist segmentation of the HE image
# b1_stat3.attrs
sopa.segmentation.stardist(b1_stat3, model_type='2D_versatile_he', image_key= 'blocco1_full_image', 
min_area=10, delete_cache=True, recover=False, prob_thresh=0.2, nms_thresh=0.6)

# let's aggregate genes counts in the bins in the new segmented nuclei polys.
sopa.aggregate(b1_stat3, key_added = 'nuclei_counts', min_transcripts=1, min_intensity_ratio=0.1)

# nuclei colorati in base all'area
plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_shapes('stardist_boundaries', color = 'area', fill_alpha = 1, outline_alpha = 1).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_nuclei_mask.png')

# we will have 2 new tables, one with nuclei_counts and one with mean intensity(segmentation image) per nuclei

