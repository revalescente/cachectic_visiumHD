import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import pandas as pd
import sopa
import json
import os
from shapely.affinity import scale
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import py_scripts.pp_sdata.pp_functions as pp

zarr_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples"
#zarr_path = "/mnt/europa/valerio/data/zarr_store/samples"
sdata = sd.read_zarr(f"{zarr_path}/blocco1_sham")

# --- Configuration ---
# Adjust these keys based on your specific sdata object structure

# SPACERANGER
COORD_SYS = 'blocco1'           # From your sdata printout
TABLE_KEY = 'segmentation_counts'
IMAGE_KEY = f'{sdata.path.stem}_full_image'        # The background image
SHAPES_KEY = f'{sdata.path.stem}_nuclei_boundaries'  # The shapes matching your table (28506 count match)
#My Segmentation
# COORD_SYS = 'blocco1'           # From your sdata printout
# TABLE_KEY = 'nuclei_counts'
# IMAGE_KEY = 'full_image'        # The background image
# SHAPES_KEY = 'filtered_nuclei'  # The shapes matching your table (28506 count match)
#' I need to sample single nucleus from specific cell types and plot them to analyze possible segmentation errors
#' I need to sample the cell_id of the nuclei for the various cell_types, get the extent, bounding box and then plot the extracted shapes 

# nuclei filtering
print("This is a spaceranger dataset so we need to filter the nucleus area %in% [5, 100]")
px_to_um = 0.316
area_conversion = px_to_um ** 2  # = 0.099856
min_area = 5
max_area = 100
print(f"Starting with {len(sdata[SHAPES_KEY])} shapes and {sdata[TABLE_KEY].n_obs} cells in table.")
print("Calculating areas...")
#sdata[SHAPES_KEY] = sdata[SHAPES_KEY].set_index('cell_id')
nuclei_pols = sdata.transform_element_to_coordinate_system(SHAPES_KEY, target_coordinate_system = COORD_SYS)
area_um = nuclei_pols.geometry.area * area_conversion
# 2. Add area to the AnnData table
# We use reindex to ensure we match the correct area to the correct cell_id in the table
sdata[TABLE_KEY].obs['area_um'] = area_um.reindex(sdata[TABLE_KEY].obs.index)
# 3. Create a Filter Mask
# Keep area >= 5 AND area <= 100
# (This removes area < 5 OR area > 100)
mask = (sdata[TABLE_KEY].obs['area_um'] >= min_area) & (sdata[TABLE_KEY].obs['area_um'] <= max_area)
# Check how many we are keeping
n_before = sdata[TABLE_KEY].n_obs
n_after = mask.sum()
print(f"Filtering: Keeping {n_after} / {n_before} nuclei ({n_before - n_after} removed)")
# 4. Apply Filter to the Table (AnnData)
sdata[TABLE_KEY] = sdata[TABLE_KEY][mask].copy()
# 5. Apply Filter to the Shapes (GeoDataFrame)
# We subset the shapes to match the IDs remaining in the filtered table
valid_ids = sdata[TABLE_KEY].obs.index
sdata[SHAPES_KEY] = sdata[SHAPES_KEY].loc[valid_ids]



print("Filtering complete.")

# Extraction of cell_id --------------------------------------------------------

# 1. Define the target cell types
which_cells = ['MuSC', 'Immune_Cells', 'Myonuclei_NMJ', 'Myonuclei_IIb']

# 2. Access the obs DataFrame
obs_df = sdata[TABLE_KEY].obs

# 3. Initialize a list to store the sampled cell_ids
collected_samples = []
# 4. Loop and sample
for cell_type in which_cells:
    # Filter for the specific cell type
    subset = obs_df[obs_df['labels_singler'] == cell_type]
    n_available = len(subset)
    if n_available > 0:
        # Sample up to 5
        n_sample = min(5, n_available)
        sampled_subset = subset.sample(n=n_sample, random_state=42)
        # Append the actual rows to our list
        collected_samples.append(sampled_subset)
        print(f"Sampled {n_sample} cells for {cell_type}")


# 5. Concatenate into a single DataFrame
summary_df = pd.concat(collected_samples)

# Plotting of the bounding box of the nuclei
sdata[SHAPES_KEY].set_index('cell_id', inplace=True)

#OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/my_segmentation/prova_annotazione'
OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/spaceranger/cell_types/b1_sham_buffer=100'
BUFFER = 100 

# --- Loop through the sampled cells ---
for index, row in summary_df.iterrows():
    cell_id = row['cell_id']
    cell_type = row['labels_singler']
    save_path = os.path.join(OUTPUT_DIR, f"{cell_type}_{cell_id}.png")
    print(f"Processing {cell_id} ({cell_type})...")
    try:
        # 1. Isolate the specific shape for this cell
        # We create a temporary element in sdata to render just this one cell in red
        temp_shape_name = "temp_highlight_shape"
        sdata[temp_shape_name] = sdata[SHAPES_KEY].loc[[cell_id]]
        # 2. Calculate the extent (bounding box) of this single cell
        extent = sd.get_extent(sdata[temp_shape_name], coordinate_system=COORD_SYS)
        # Extract coordinates
        min_x, max_x = extent['x']
        min_y, max_y = extent['y']
        # 3. Define the plot window (Cell bbox + buffer)
        bbox_min = [min_x - BUFFER, min_y - BUFFER]
        bbox_max = [max_x + BUFFER, max_y + BUFFER]
        # 4. Setup Plot
        plt.figure(figsize=(10, 10))
        ax = plt.gca()
        # 5. Perform the Query and Render
        # We crop the sdata to the bounding box first for efficiency
        sdata.query.bounding_box(
            axes=["x", "y"],
            min_coordinate=bbox_min,
            max_coordinate=bbox_max,
            target_coordinate_system=COORD_SYS,
        ).pl.render_images(
            IMAGE_KEY
        ).pl.render_shapes(
            SHAPES_KEY, 
            outline=True, 
            outline_alpha=0.5, 
            outline_width=1.0, 
            fill_alpha=0
        ).pl.render_shapes(
            temp_shape_name, 
            fill_alpha=0, 
            outline=True, 
            outline_color="yellow",
            outline_width=2,
            outline_alpha=1
        ).pl.show(
            ax=ax, 
            title=f"{cell_type}\n{cell_id}", 
            coordinate_systems=COORD_SYS,
            save = save_path
        )
        plt.close()
        # Cleanup the temporary shape from sdata
        del sdata[temp_shape_name]
    except KeyError:
        print(f"  Error: Cell ID {cell_id} not found in shape '{SHAPES_KEY}'.")
    except Exception as e:
        print(f"  Error plotting {cell_id}: {e}")

# ---------------------------------------------------------------------

# extract the cell_id of the nuclei to plot:
cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']

# A list to hold all the results
all_extremes = []


# Loop through each feature and get the extremes
for column in cols:
    if column in sdata['nuclei_counts'].obs.columns:
        top_3, bottom_3 = get_extreme_observations(sdata['nuclei_counts'].obs, column, number_observations=3)
        
        # Add a column to identify the feature and extreme type
        top_3 = top_3.assign(feature=column, extreme_type='Top 3')
        bottom_3 = bottom_3.assign(feature=column, extreme_type='Bottom 3')
        
        all_extremes.extend([top_3, bottom_3])

# Combine everything into a single DataFrame for a nice summary
summary_df = pd.concat(all_extremes)

# Reorder columns for better readability
display_cols = ['feature', 'extreme_type', 'cell_id'] + feature_columns
summary_df = summary_df.reset_index().rename(columns={'index': 'cell_id'})
summary_df = summary_df[display_cols]


# Display the final summary table
print(summary_df.to_string())

# we have the info about the cell_id of the extremes nuclei. 
# now let's try to plot one nucleus

# 1. Define the list of cell_ids you want to keep.
#    Replace these example IDs with your actual list.
ids_to_keep = summary_df.index

# 2. Filter the GeoDataFrame using .loc
#    This selects only the rows whose index (cell_id) is in your list.
nucleo_geo = sdata['blocco4_nuclei_boundaries'].loc[[summary_df.index[0]]]
sdata['nucleo_shape'] = nucleo_geo

sd.get_extent(sdata['nucleo_shape'], coordinate_system='blocco4', exact=True)
# {'x': (3955.5047917988886, 3966.4837517193064), 'y': (12629.535219786241, 12649.252859600709)}

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[3955.5047917988886-50, 12629.535219786241-50],
    max_coordinate=[3966.4837517193064+50, 12649.252859600709+50],
    target_coordinate_system="blocco4",
).pl.render_images("blocco4_full_image").pl.render_shapes("blocco4_nuclei_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.render_shapes("nucleo_shape", color = "red", fill_alpha = 0.5).pl.show(ax = ax, coordinate_systems="blocco4", save = 'output_python/nuclei_Strange/most_ecc_1.png')


# nuclei area histograms ----

# Conversion factor
px_to_um = 0.316
area_conversion = px_to_um ** 2  # = 0.099856

# Convert areas to μm²

area_micron = calculated_areas * area_conversion

# min_area=5/area_conversion
# max_area=100/area_conversion
# 2. Create the plot (do not use plt.show())
plt.figure(figsize=(10, 6))
plt.hist(area_micron, bins=50, color='skyblue', edgecolor='black', alpha=0.7)

# Add 2 vertical lines
# plt.axvline(min_area, color='red', linestyle='--', linewidth=2, label=f"min_area = {min_area}")
# plt.axvline(max_area, color='orange', linestyle='--', linewidth=2, label=f"max_area = {max_area}")

# Add labels and legend
plt.title(f"Distribution of Nuclei Area (n={len(area_micron)})")
plt.xlabel("Area in micron")
plt.ylabel("Frequency")
plt.grid(axis='y', alpha=0.5)
plt.legend()

# 3. Save to disk
output_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/nuclei_area_histogram_sr_b1_sham.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')

# 4. Close the figure to free memory and prevent display
plt.close()

print(f"Histogram saved to {output_path}")
