import os
import spatialdata as sd
import geopandas as gpd
from pathlib import Path
from py_scripts.utils.utils_fun import read_from_json

# --- Configuration ---
# 1. Info about the samples
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

# 2. Directory of the sdata
spatialdata_dir = "/mnt/europa/valerio/data/zarr_store/samples"

# 3. Directory to save the extracted GeoDataFrames
output_dir = "/mnt/europa/valerio/data/extracted_shapes"
os.makedirs(output_dir, exist_ok=True)

# 4. The specific key of the shape element you want to extract
# Check sdata.shapes.keys() if you are unsure (e.g., 'nuclei', 'cell_boundaries', 'circles')
shape_key = "nuclei_boundaries" 

# --- Processing Loop ---
for blocco, samples in samples_dict.items():
    for sample, _ in samples.items():
        sample_id = f"{blocco}_{sample}"
        spatialdata_path = os.path.join(spatialdata_dir, f"{sample_id}.zarr")
        print(f"Processing {sample_id}...")
        try:
            # 1. Read the SpatialData object
            # We can use mode='r' since we are only reading
            sdata = sd.read_zarr(spatialdata_path)
            # 2. Check if the shape key exists
            if shape_key not in sdata.shapes:
                print(f"  -> Warning: Shape key '{shape_key}' not found in {sample_id}")
                continue
            shapes_dict, _ = sd.match_element_to_table(sdata, shape_key, "nuclei_counts")
            gdf = shapes_dict[shape_key]
            # 3. Extract the GeoDataFrame
            # This returns a geopandas.GeoDataFrame
            gdf['name'] = gdf.index
            # 4. (Optional) Cleaning
            # Sometimes SpatialData keeps the index named 'instance_id' or similar.
            # Resetting index can make it safer for some file formats if the index is complex.
            # gdf = gdf.reset_index()
            # 5. Save to disk
            # Parquet is the best format for GeoDataFrames (fast, preserves types)
            output_path = os.path.join(output_dir, f"{sample_id}_{shape_key}.geojson")
            gdf.to_file(output_path, driver="GeoJSON")
            print(f"  -> Extracted {len(gdf)} shapes to {output_path}")
        except Exception as e:
            print(f"  -> Error processing {sample_id}: {e}")



print("Done.")
