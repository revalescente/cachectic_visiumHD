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
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf
#
from importlib import reload
# reload(sf)


# sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26foxO.zarr")
# sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26.zarr")
sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26SMAD23.zarr")

# control the attributes
# sopa.utils.set_sopa_attrs(
#     sdata,
#     cell_segmentation_key="blocco4_full_image",
#     tissue_segmentation_key="blocco4_full_image",
#     bins_table_key="filtered"
# )

# To limit the segmentation to the area of the image within the tissue we have to give the "intissue" shape 
# the name "region_of_interest". In the next version it's possibile will be added the possibility to set 
# the "intissue" shape's name as an attrs.

# define the region of interest aka in tissue area - must be this name to be seen in the segmentation function
sdata['region_of_interest'] = sdata['blocco4_intissue'][sdata['blocco4_intissue'].name=="c26SMAD23"]

sopa.make_image_patches(sdata, patch_width = 500)

sopa.segmentation.stardist(sdata, model_type='2D_versatile_he', 
image_key= sdata.attrs['cell_segmentation_image'], min_area=10, delete_cache=True, 
recover=False, prob_thresh=0.2, nms_thresh=0.6, key_added = 'blocco4_nuclei_boundaries')

# It worked, we only have segmentation results from inside the tissue
plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images('blocco4_full_image', scale = "scale2").pl.render_shapes("blocco4_intissue", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.render_shapes("blocco4_nuclei_boundaries", color = "green", outline=True, outline_alpha=0, fill_alpha=1
).pl.show(ax = ax, coordinate_systems="blocco4", save = 'output_python/roba_da_buttare/sdata_nuclei_b4_c26smad23.png')

# ------------------------------------------------------------------------------

# Before aggregation we need to filter all low quality nuclei
# We use morphological characteristics of the nuclei to that:
# area
# eccentricity: aka how much spheric is the nucleus? if ecc = 0 then it's a sphere
#               if ecc = 1 it's a parabola 
# solidity: Is it convex or concave? if sol = 1 it's convex,
#           concave otherwise
# extent: how much squared is the nucleus? ratio between nucleus area and bbox area

#
# Features extraction
try:
    # Assuming 'sdata' is your loaded SpatialData object
    features_df = sf.features_extraction(sdata, nuclei_element_name="blocco4_nuclei_boundaries")
    print("Successfully extracted features:")
    print(features_df.head())
except ValueError as e:
    print(f"Error: {e}")

cols = ['area', 'eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']

for col in cols:
    plt.figure(figsize=(20, 20))
    plt.hist(features_df[col], bins=50, color='skyblue', edgecolor='black')
    plt.title(f'Histogram of {col}')
    plt.xlabel(col)
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f'figures/output_python/{col}_histogram.png')  # Saves the plot as PNG
    plt.close()  # Close the figure to avoid display overlap

# -------
# let's aggregate genes counts in the bins in the new segmented nuclei polys.
sopa.aggregate(sdata, key_added = 'nuclei_counts', bins_key= "filtered", 
shapes_key = "blocco4_nuclei_boundaries", expand_radius_ratio=1, min_transcripts=1, 
min_intensity_ratio=0.1)


# add values inside the .obs of the corresponding table
cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
sdata['nuclei_counts'].obs = sdata['nuclei_counts'].obs.join(features_df[cols], how='left')

# plotting different nuclei ----------------------

# collect info about the most interesting nuclei
cols_to_analyze = ['area', 'eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']

# This will now produce a DataFrame with no duplicate columns
extreme_nuclei_df = sf.find_extreme_observations(
    sdata=sdata,
    table_key='nuclei_counts',
    feature_columns=cols_to_analyze,
    n_extremes=3
)

# Verify the columns are clean
print(extreme_nuclei_df.columns)


# Visualize these nuclei

ids_to_visualize = extreme_nuclei_df['cell_id']
# 
# Call the function to generate and save the plots
sf.plot_nuclei(
    sdata = sdata,
    extreme_nuclei_df = extreme_nuclei_df,
    shapes_key = "blocco4_nuclei_boundaries",
    save_dir = 'output_python/nuclei_Strange/',
    padding = 50,
    figsize = (15, 15))


# Filtering nuclei -------------------------------------------------------------
filters = {'area': (10, 4000),  # eliminiamo i più piccoli e quelli troppo grandi
'eccentricity': (None, 0.9), # eliminiamo i nuclei più parabolici (di solito con un lato simil retto)
'solidity': (0.7, None), # eliminiamo i più concavi
'extent': (0.2, None)} # eliminiamo i nuclei più irregolari (un misto tra concavi e parabolici)

filtered = sf.nuclei_filtering(features_df, filters)

# Get the cell_ids to keep
cell_ids_to_keep = filtered.index.tolist()

# Filter the GeoDataFrame by index
sdata['blocco4_filtered_nuclei'] = sdata['blocco4_nuclei_boundaries'].loc[cell_ids_to_keep]

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images("blocco4_full_image").pl.render_shapes("blocco4_filtered_nuclei", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco4", save = 'output_python/b4_c26_filtered_nuclei.png')


# ------------------------------------------------------------------------------

# assign bins to nuclei
sdata['blocco4_bins_nuclei'] = sf.assign_bins_to_nearest_nucleus(
sdata['blocco4_intissue_002um'], sdata['blocco4_filtered_nuclei'])


# Plots:
# intensities per channel per nucleus.
# Assuming sdata['nuclei_counts'].obsm['intensities'] is a pandas DataFrame
df = sdata['nuclei_counts'].obsm['intensities']

fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns

channels = ['r', 'g', 'b']
colors = ['red', 'green', 'blue']
titles = ['Red channel', 'Green channel', 'Blue channel']

for ax, channel, color, title in zip(axs, channels, colors, titles):
    ax.hist(df[channel], bins=30, color=color, alpha=0.8)
    ax.set_title(title)
    ax.set_xlabel('Intensity')
    ax.set_ylabel('Count')

plt.tight_layout()
plt.savefig('/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/roba_da_buttare/nuclei_rgb_intensities.png')
plt.close()

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images('blocco4_full_image', scale = "scale2").pl.render_shapes("blocco4_intissue", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.render_shapes("blocco4_nuclei_boundaries", color = "green", outline=True, outline_alpha=0, fill_alpha=1
).pl.show(ax = ax, coordinate_systems="blocco4", save = 'output_python/roba_da_buttare/sdata_nuclei_b4.png')

# ------------------------------------------------------------------------------
# grafici di controllo: rasterize
element_extent = sd.get_extent(sdata['blocco9_full_image'], coordinate_system='blocco9', exact=False)

sdata['blocco9_rast_nuclei'] = sd.rasterize(
  sdata['blocco9_nuclei_boundaries'],
  ["x", "y"],
  min_coordinate=[element_extent['x'][0],element_extent['y'][0]],
  max_coordinate=[element_extent['x'][1],element_extent['y'][1]],
  target_coordinate_system='blocco9',
  target_unit_to_pixels=1,
)


# plotting the rasterized data
plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[6000, 4000],
    max_coordinate=[8000, 5000],
    target_coordinate_system="blocco9",
)sdata.pl.render_images("blocco9_rast_nuclei").pl.render_shapes("blocco9_nuclei_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco9", save = 'output_python/roba_da_buttare/sdata_nuclei.png')


