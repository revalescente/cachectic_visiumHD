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

# let's display the areas where no expression is detected as transparent
new_cmap = set_zero_in_cmap_to_transparent(cmap="viridis")
new_cmap

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
sopa.aggregate(b1_stat3, key_added = 'nuclei_counts', expand_radius_ratio=1, min_transcripts=1, min_intensity_ratio=0.1)

# nuclei colorati in base all'area
plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_shapes('stardist_boundaries', color = 'area', fill_alpha = 1, outline_alpha = 1).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_nuclei_mask.png')

# we will have 2 new tables, one with nuclei_counts and one with mean intensity(segmentation image) per nuclei

# I have now a new table "nuclei_counts" result of the segmentation result, aggregation of the counts in the bins but 
# somehow we lost the annotation between bins and cells, which bins belong to the cells? why there isn't an index between them??

# let's redo the aggregation step and specify everything meaningful
sopa.aggregate(b1_stat3, key_added = 'nuclei_counts', bins_key= "cleaned", shapes_key = "stardist_boundaries", expand_radius_ratio=1, min_transcripts=1, min_intensity_ratio=0.1)

# nothing meaninful is different

# qc info
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
b1_stat3['nuclei_counts'].var["mt"] = b1_stat3['nuclei_counts'].var_names.str.startswith("Mt-")
# ribosomal genes
b1_stat3['nuclei_counts'].var["ribo"] = b1_stat3['nuclei_counts'].var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
b1_stat3['nuclei_counts'].var["hb"] = b1_stat3['nuclei_counts'].var_names.str.contains("^Hb[^(p)]")
sc.pp.calculate_qc_metrics(
    b1_stat3['nuclei_counts'], qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# some plots 
# area < 1000
adata = b1_stat3['nuclei_counts'][b1_stat3['nuclei_counts'].obs['area'] < 1000, : ]


# Get numeric columns only
numeric_columns = adata.obs.select_dtypes(include='number').columns.tolist()

# Split columns into two roughly equal groups
split = len(numeric_columns) // 2
group1 = numeric_columns[:split]
group2 = numeric_columns[split:]

def plot_grid(columns, filename, ncols=4):
    nrows = int(np.ceil(len(columns) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
    axes = axes.flatten()
    for i, col in enumerate(columns):
        ax = axes[i]
        adata.obs[col].hist(ax=ax, bins=50)
        ax.set_title(col)
    # Hide unused subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close(fig)

# Create and save the two grid plots
plot_grid(group1, "b1_stat3_obs_hist_group1.png", ncols=4)
plot_grid(group2, "b1_stat3_obs_hist_group2.png", ncols=4)

# end segmentation
#------------------------------------------------------------------------------------------------------------------
# start cleaning segmentation and understanding of what's what and where.

# let's try to annotate the nuclei_counts with the bins, so let's match nuclei and bins

sd.match_element_to_table(b1_stat3, "stardist_boundaries", "nuclei_counts") # not possible because annotation is lost with aggregation

# let's try to filter element with element "intissue_002um" filtered with "stardist_boundaries"
# I use the aggregate function between shapes
# NO! it aggregate some finite value (like gene expression (1 gene))

# let's try sjoin of sopa
# i need to now which bin belong to which nuclei
# Want to merge to get each bin's nuclei. 
b1_stat3['intissue_002um'].head() 
b1_stat3['stardist_boundaries'].head() # the important thing is the index column, maybe I should move it in a column
b1_stat3['stardist_boundaries']['cell_id']

# bins_with_nuclei = sopa.spatial.sjoin(b1_stat3, "stardist_boundaries", "intissue_002um", how='inner', target_coordinate_system="blocco1")


# # 
# # extract bins shapes keeping the index 'location_id' for the filtering
# bins = b1_stat3['intissue_002um']  # location_id becomes a column
# bins['location_id'] = bins['location_id'].astype('categorical') # trasformiamo in categoriale
# nuclei_bounds = b1_stat3['stardist_boundaries']
# nuclei_bounds['cell_id'] = nuclei_bounds.index
# 
# # let's do as copilot say: filtering bins inside the cells mantaining the cell_id info
# filtered_bins = gpd.sjoin(
#     bins,
#     nuclei_bounds[['geometry', 'cell_id']],
#     how='inner',
#     predicate='intersects'  # or 'within'
# )

# Extract bins' cell_id informations

# extract bins shapes keeping the index 'location_id' for the filtering
bins = b1_stat3['intissue_002um']  # location_id becomes a column
bins['location_id'] = bins['location_id'].astype('category') # trasformiamo in categoriale
nuclei_bounds = b1_stat3['stardist_boundaries']
nuclei_bounds['cell_id'] = nuclei_bounds.index

# Spatial join to the all intersecting bin-nucleus pairs, mantaining the cell_id info 
filtered_bins = gpd.sjoin(
    bins,
    nuclei_bounds[['geometry', 'cell_id']],
    how='inner',
    predicate='intersects'  # or 'within'
)

# add bounds geometries to the filtered_bins - cell_id is a column in nuclei_bounds:
filtered_bins = filtered_bins.merge(
    nuclei_bounds[['cell_id', 'geometry']].rename(columns={'geometry': 'nucleus_geometry'}),
    on='cell_id'
)
# assign the bins to the nearest nucleus boundary ------------------------------

# compute bin centroid
filtered_bins['bin_centroid'] = filtered_bins.geometry.centroid
# add distance column between centroids and borders 
filtered_bins['distance'] = filtered_bins.apply(
    lambda row: row['bin_centroid'].distance(row['nucleus_geometry']),
    axis=1
)

# for each bin, keep only the closest nucleus
# This will keep only one row per bin (location_id), corresponding to the nearest intersected nucleus
nearest_bins = filtered_bins.loc[
    filtered_bins.groupby('location_id', observed=True)['distance'].idxmin()
].reset_index(drop=True)

# Step 5: Optional - select relevant columns for output
# For example, keep bin id, assigned nucleus id, and geometries
bins_nuclei = nearest_bins[['location_id', 'cell_id', 'name', 'geometry', 'distance']]

# result is now a GeoDataFrame with each bin assigned to its closest intersecting nucleus boundary


# add in the spe object
intersection_parse = ShapesModel.parse(filtered_bins), transformations={'blocco1': Identity()})
b1_stat3['bins_nuclei2'] = intersection_parse
ShapesModel().validate(intersection_parse)

# # fix the cell_id of the anndata with the cell_id of the bins
# cell_ids_to_keep = b1_stat3['bins_nuclei']['cell_id'].unique()
# b1_stat3['nuclei_counts'] = b1_stat3['nuclei_counts'][b1_stat3['nuclei_counts'].obs['cell_id'].isin(cell_ids_to_keep), :]
# # add the location_id of the bins_nuclei to annotate the shapes
# 
# # Get cell_id-location_id mapping (one location_id per cell_id)
# location_map = (
#     b1_stat3['bins_nuclei'][['cell_id', 'location_id']]
#     .drop_duplicates('cell_id') # if there are duplicates, take the first
#     .set_index('cell_id')['location_id']
# )
# 
# # Assign new column in filtered.obs
# filtered.obs['location_id'] = filtered.obs['cell_id'].map(location_map)
# 
# print(f"Number of missing location_id: {filtered.obs['location_id'].isna().sum()}")
# 
# b1_stat3['nuclei_countsB']  = TableModel.parse(filtered)

# fare grafici per vedere che abbiamo in mano

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_images("blocco1_conv_image", cmap = "grey").pl.render_shapes(
    "bins_nuclei2", outline_alpha=1, fill_alpha=1, color = "green", cmap = new_cmap
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_bins_nuclei2.png')



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


# CONTINUARE DA QUI 





# let's try to plot only the nuclei with area > 2000 and test their presence and type.

dummarea =  b1_stat3['nuclei_counts'].obs['area'] > 2000

plt.figure(figsize=(50, 50))
ax = plt.gca()
b1_stat3.pl.render_shapes('stardist_boundaries', color = dummarea , fill_alpha = 1, outline_alpha = 1).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/b1_stat3_nuclei_mask.png')
