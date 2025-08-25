import spatialdata as sd
import spatialdata_plot
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import sopa
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
from py_scripts.segmentation.segm_functions import (assign_bins_to_nearest_nucleus, nuclei_filtering) 

# let's display the areas where no expression is detected as transparent
new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
new_cmap

# read data
b1_stat3 = read_zarr_standardized(os.path.expanduser("~/data/b1_stat3.zarr"))

# set-up/check attributes
sopa.utils.set_sopa_attrs(
    b1_stat3,
    cell_segmentation_key="blocco1_full_image",
    tissue_segmentation_key="blocco1_full_image",
    bins_table_key="cleaned",
    boundaries_key = "blocco1_intissue"
)

# divide the image in patches to segments individually first (overlapping patches)
# sopa.make_image_patches(sdata, patch_width=2000, patch_overlap=50, image_key=None, key_added=None)
sopa.make_image_patches(b1_stat3, patch_width = 500, patch_overlap = 50)

# stardist segmentation of the HE image
# b1_stat3.attrs
sopa.segmentation.stardist(b1_stat3, model_type='2D_versatile_he', image_key= 'blocco1_full_image', 
min_area=10, delete_cache=True, recover=False, prob_thresh=0.2, nms_thresh=0.6)

# let's aggregate genes counts in the bins in the new segmented nuclei polys.
sopa.aggregate(b1_stat3, key_added = 'nuclei_counts', bins_key= "cleaned", 
shapes_key = "stardist_boundaries", expand_radius_ratio=1, min_transcripts=1, 
min_intensity_ratio=0.1)

#-------------------------------------------------------------------------------

# problem: 
#' we have to filter nuclei based on dimensions, shapes, ecc
#' we have to solve bins belonging to more then one nucleus:
#'        solution: assign to the nearest nucleus -            DONE!

# filter nuclei ----------------------------------------------------------------

# let's work with skimage regionprops, i need to rasterize the shapes to become labels

sd.get_extent(b1_stat3['stardist_boundaries'], coordinate_system='blocco1', exact=True)
# exemple {'x': (1994.433953335766, 15784.49036300504), 'y': (299.5143364244573, 7999.5)}

# rasterize shapes
b1_stat3["raster_nuclei"] = sd.rasterize(
    b1_stat3["stardist_boundaries"],
    ["x", "y"],
    min_coordinate=[1994.433953335766, 299.5143364244573],
    max_coordinate=[15784.49036300504, 7999.5],
    target_coordinate_system="blocco1",
    target_unit_to_pixels=1,
)
# plotting the rasterized data
plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[6000, 4000],
    max_coordinate=[8000, 5000],
    target_coordinate_system="blocco1",
).pl.render_images("raster_nuclei").pl.render_shapes("stardist_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_raster_nuclei.png')

# transform into integer values

# If your array is an xarray with shape (1, H, W), get the first channel
label_mask = b1_stat3['raster_nuclei'].values

# If shape is (1, H, W), squeeze to (H, W)
if label_mask.ndim == 3 and label_mask.shape[0] == 1:
    label_mask = label_mask[0]

# Cast to integer type
label_mask = label_mask.astype(np.int32)

# add features of the rasterized nuclei boundaries

# 2. Compute regionprops
props = regionprops_table(label_mask, properties=[
    'label', 'area', 'eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length'])
props_df = pd.DataFrame(props)

# Get mapping dictionary from xarray attributes
label_to_id = b1_stat3['raster_nuclei'].attrs['label_index_to_category']

# Map: create a new column with the nucleus ID
props_df['cell_id'] = props_df['label'].map(label_to_id)

# now we have a lot of morphological information about the nuclei shapes. 
# Let's visualize the collected info

# area: define a min - max area
# eccentricity: if = 0 then it's a sphere
# solidity: convexity measure (if = 1 it's convex)
# extent: ratio between nucleus area and bbox area
# List of columns to plot
cols = ['area', 'eccentricity', 'solidity', 'extent']

for col in cols:
    plt.figure(figsize=(20, 20))
    plt.hist(props_df[col], bins=50, color='skyblue', edgecolor='black')
    plt.title(f'Histogram of {col}')
    plt.xlabel(col)
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f'figures/output_python/{col}_histogram.png')  # Saves the plot as PNG
    plt.close()  # Close the figure to avoid display overlap


# let's filter now

filters = {'area': (10, 1000), 'eccentricity': (0.4, None), 
'solidity': (0.7, None), 'extent': (0.1,0.9)}

filtered = nuclei_filtering(props_df, filters)

# filter nuclei boundaries -----------------------------------------------------

# Get the cell_ids to keep
cell_ids_to_keep = filtered['cell_id'].tolist()

# Filter the GeoDataFrame by index
prova = b1_stat3['stardist_boundaries'].loc[cell_ids_to_keep]

b1_stat3['nuclei_boundaries'] = prova

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images("raster_nuclei").pl.render_shapes("nuclei_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_filtered_nuclei.png')


# plotting the rasterized data
plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[6000, 4000],
    max_coordinate=[8000, 5000],
    target_coordinate_system="blocco1",
).pl.render_images("raster_nuclei").pl.render_shapes("nuclei_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_filtered_nuclei_zoom.png')

# ------------------------------------------------------------------------------

# Assign bins to nuclei

# Example usage:
# sdata['bin_shapes'] = assign_bins_to_nearest_nucleus(b1_stat3['intissue_002um'], b1_stat3['stardist_boundaries'])
# 
b1_stat3['bins_nuclei'] = assign_bins_to_nearest_nucleus(b1_stat3['intissue_002um'], b1_stat3['stardist_boundaries'])

# b1_stat3['bins_nuclei2']['cell_id'] = b1_stat3['bins_nuclei2']['cell_id'].astype('category')

# plot of the resulting bins - with and without zooming
plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "bins_nuclei2", outline_alpha=1, fill_alpha=1, color = "green", cmap = new_cmap
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_bins_nuclei.png')



plt.figure(figsize=(50, 50))
ax = plt.gca()

b1_stat3.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[6000, 4000],
    max_coordinate=[8000, 5000],
    target_coordinate_system="blocco1",
).pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "bins_nuclei2", outline_alpha=1, outline_width=0.5, fill_alpha=1, color = "green", cmap = new_cmap
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_bins_nuclei2.png')



# let's make a plot that has:
# blocco1_conv_image, bins_nuclei colored by cell_id, nuclei boundaries not filled 

plt.figure(figsize=(50, 50))
ax = plt.gca()

b1_stat3.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[6000, 4000],
    max_coordinate=[8000, 5000],
    target_coordinate_system="blocco1",
).pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "bins_nuclei2", outline_alpha=0.5, outline_width=0.5, fill_alpha=1, color = "green", cmap = new_cmap
).pl.render_shapes(
    "stardist_boundaries", outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_bins_nuclei_both.png')

