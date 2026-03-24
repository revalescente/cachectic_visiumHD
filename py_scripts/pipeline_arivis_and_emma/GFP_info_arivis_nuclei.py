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

# Automatically add repo root to path to fix "No module named 'py_scripts'"
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
if repo_root not in sys.path:
    sys.path.append(repo_root)

from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp
import py_scripts.segmentation.segm_functions as sf

# 1. Get the blocco key from the command line
if len(sys.argv) < 2:
    print("Usage: python process_single_blocco.py <BLOCCO_KEY>")
    sys.exit(1)

BLOCCO_KEY = sys.argv[1]

# disable autosave of sopa
sopa.settings.auto_save_on_disk = False

# paths of interest
sdata_dir = "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins"

# Load dictionary
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

if BLOCCO_KEY not in samples_dict:
    print(f"[ERROR] {BLOCCO_KEY} not found in dictionary.")
    sys.exit(1)

block_samples = samples_dict[BLOCCO_KEY]

print(f"\n{'='*50}")
print(f"=== PROCESSING BLOCK: {BLOCCO_KEY} ===")
print(f"{'='*50}")

for sample_name, sample_info in block_samples.items():
    print(f"\n--- Processing sample: {BLOCCO_KEY} - {sample_name} ---")
    sdata = sd.read_zarr(f"{sdata_dir}/{BLOCCO_KEY}_{sample_name}")
    
    polys_key = f"intissue_GFP_poly_{sample_name}"
    nuclei_key = "nuclei_arivis_poly"
    table_name = "arivis_nuclei_table"
    
    # 1. Transform elements
    left_element = sdata.transform_element_to_coordinate_system(nuclei_key, BLOCCO_KEY)
    right_element = sdata.transform_element_to_coordinate_system(polys_key, BLOCCO_KEY)
    
    # 2. FILTER: Keep only 'fibre_trattate' and 'infiammazione' (exclude sample_name)
    right_element_filtered = right_element[
        right_element['name'].str.contains('fibre_trattate|infiammazione', na=False, case=False)
    ]
    # Check if there are any polygons left after filtering
    if right_element_filtered.empty:
        print(f"    [SKIP] No 'fibre_trattate' or 'infiammazione' polygons found for {sample_name}. Moving to next sample.")
        continue
    
    # 3. Spatial Join
    joined_bins = gpd.sjoin(left_element, right_element_filtered, how='inner', predicate='intersects')
    
    table = sdata.tables[table_name].copy()
    instance_key = table.uns['spatialdata_attrs']['instance_key']
        
    # Initialize columns
    for col in ['in_treatment', 'to_discard']:
        table.obs[col] = False
            
    # Assign True to matched nuclei
    nuclei_per_poly = joined_bins.groupby('name').apply(lambda x: x.index.tolist())
    for source_key in nuclei_per_poly.index:
        col_name = None
        if "fibre_trattate" in source_key: col_name = "in_treatment"
        elif "infiammazione" in source_key: col_name = "to_discard"
            
        if col_name:
            nuclei_list = nuclei_per_poly[source_key]
            mask = table.obs[instance_key].isin(nuclei_list)
            if mask.any():
                table.obs.loc[mask, col_name] = True
                
    # 4. Calculate GFP values (Fixed variables to match current script)
    print("    Calculating GFP values for Arivis nuclei...")
    channel_aggregation = sopa.aggregation.aggregate_channels(
        sdata,                                 # changed from sdata_bbox
        image_key=f'{BLOCCO_KEY}_fluo_image', 
        shapes_key=nuclei_key,                 # changed from shape_key
        expand_radius_ratio=0, 
        mode='max', 
        no_overlap=False
    ) 
    max_values_vector = channel_aggregation.max(axis=1)
    table.obs['GFP_value'] = max_values_vector
    
    # 5. Save the updated table back to the Zarr store
    sdata.tables[table_name] = table
    # write_element securely overwrites just the table inside your existing Zarr directory
    sdata.delete_element_from_disk(table_name)
    sdata.write_element(table_name)
    
    print(f"    Finished updating {table_name} for {sample_name} and saved to disk.")
