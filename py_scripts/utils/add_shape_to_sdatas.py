import sys
import os
import argparse
import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity, get_transformation
from skimage.measure import regionprops_table
import sopa

# --- PATH SETUP ---
# Ensure py_scripts can be imported if running from within the folder
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) 
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

try:
    from utils.utils_fun import read_from_json
except ImportError:
    print("Warning: Could not import 'read_from_json'. Make sure environment is correct.")
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Add shape to each sdata's sample.")
    
    # Defaults
    default_json = "/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json"
    default_geojson = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
    default_sdata = "/mnt/europa/valerio/data/zarr_store/samples"
    
    parser.add_argument("--json", default=default_json, help="Path to sample JSON")
    parser.add_argument("--geojson_dir", default=default_geojson, help="Path to GeoJSON dir")
    parser.add_argument("--sdata_dir", default=default_sdata, help="Path to sdata dir")
    parser.add_argument("--shape_name", default="GFP_poly", help="Name of the new shape element")
    
    args = parser.parse_args()

    # 1. Define the protected names list
    protected_names = [
        'filtered_bins', 
        'filtered_nuclei', 
        'intissue_002um', 
        'intissue_poly', 
        'nuclei_boundaries'
    ]

    # 2. Check if the chosen name is protected
    if args.shape_name in protected_names:
        print(f"  [ERROR] '{args.shape_name}' is a protected core element.")
        print(f"  Please choose a different --shape_name to avoid overwriting original data.")
        return # Stop execution
      
    samples_dict = read_from_json(args.json)
    
    for block_id, conditions in samples_dict.items():
        for condition_id, details in conditions.items():
            sample_key = details['sample_key']
            print(f"\nProcessing sample: {sample_key}")
            
            # 1. Load Sdata
            zarr_path = os.path.join(args.sdata_dir, f"{sample_key}.zarr")
            if not os.path.exists(zarr_path):
                print(f"  [SKIP] Zarr not found: {zarr_path}")
                continue
            
            sdata = sd.read_zarr(zarr_path)
            
            if args.shape_name in sdata.shapes:
                print(f"  [Alert] '{args.shape_name}' already exists. Deleting from disk and overwriting in memory...")
                # Use the spatialdata utility to clear the Zarr store entries
                sdata.delete_element_from_disk(args.shape_name)

            
            # 2. Get the transformation of the reference full image
            transform = get_transformation(sdata["full_image"], get_all = True)

            # 3. Load GeoJSON
            geojson_file = os.path.join(args.geojson_dir, f"{sample_key}.geojson")
            polys_loaded = False
            
            if os.path.exists(geojson_file):
                try:
                    sample_poly = gpd.read_file(geojson_file)
                    if not sample_poly.empty:
                        sample_poly = sample_poly.set_crs(None, allow_override=True)
                        sample_poly_parsed = ShapesModel.parse(sample_poly, transformations=transform)
                        sdata[args.shape_name] = sample_poly_parsed
                        polys_loaded = True
                    else:
                        print(f"  [WARN] GeoJSON empty: {sample_key}")
                except Exception as e:
                    print(f"  [ERROR] GeoJSON error: {e}")
            else:
                print(f"  [INFO] No GeoJSON found for {sample_key}")

            # 8. Save
            sdata.write_element(args.shape_name)
            print(f"  [SUCCESS] Saved '{args.shape_name}'")

if __name__ == "__main__":
    main()
