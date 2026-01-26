import sys
import os
import argparse
import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity
from skimage.measure import regionprops_table
import sopa

# --- PATH SETUP ---
# Ensure py_scripts can be imported if running from within the folder
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir) 
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

try:
    from py_scripts.utils.utils_fun import read_from_json
except ImportError:
    print("Warning: Could not import 'read_from_json'. Make sure environment is correct.")
    sys.exit(1)


def compute_morphology(sdata, nuclei_key, coord_sys):
    """
    Computes morphological features by rasterizing the nuclei shapes.
    """
    print(f"    Computing morphological features for {nuclei_key}...")
    raster_key = f"raster_{nuclei_key}"
    
    # 1. Get extent
    try:
        element_extent = sd.get_extent(sdata[nuclei_key], coordinate_system=coord_sys, exact=True)
    except Exception as e:
        print(f"      [WARN] Could not get extent: {e}. Skipping morphology.")
        return None

    # 2. Rasterize
    # Note: For Shapes, this creates a raster. If shapes overlap, usually last wins.
    # We hope sd.rasterize preserves indices in the values or we rely on the implementation details.
    try:
        sdata[raster_key] = sd.rasterize(
            sdata[nuclei_key],
            axes=["x", "y"],
            min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
            max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
            target_coordinate_system=coord_sys,
            target_unit_to_pixels=1,
        )
    except Exception as e:
        print(f"      [WARN] Rasterization failed: {e}. Skipping morphology.")
        return None

    # 3. Prepare Mask
    # Squeeze to (H, W) and ensure int
    try:
        label_mask = sdata[raster_key].data.squeeze().astype(np.int32)
        # Handle dask array if necessary
        if hasattr(label_mask, 'compute'):
            label_mask = label_mask.compute()
    except Exception as e:
        print(f"      [WARN] Mask preparation failed: {e}")
        return None

    # 4. Compute Features (regionprops)
    properties_to_extract = [
        'label', 'area', 'eccentricity', 'solidity', 'extent', 
        'major_axis_length', 'minor_axis_length'
    ]
    try:
        props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))
    except Exception as e:
        print(f"      [WARN] regionprops failed (mask might be empty): {e}")
        return None

    if props_df.empty:
        print("      [WARN] Morphology DataFrame is empty.")
        return None

    # 5. Map Labels to Cell IDs
    # Strategy A: Check for explicit mapping in attrs (User's original method)
    mapped = False
    if 'label_index_to_category' in sdata[raster_key].attrs:
        try:
            label_to_id = sdata[raster_key].attrs['label_index_to_category']
            props_df['cell_id'] = props_df['label'].map(label_to_id)
            mapped = True
        except Exception:
            pass
    
    # Strategy B: If no explicit mapping, assume label integer k corresponds to shape at index k-1 
    # (or k if 0-based, but regionprops usually starts at 1)
    if not mapped:
        # Get the actual index from the shapes dataframe
        try:
            nuclei_index = sdata[nuclei_key].index
            # Create a mapper: label 1 -> index[0], label 2 -> index[1]...
            # This assumes rasterization used standard encoding (1-based index)
            # We filter props to only labels that exist within the range of the dataframe
            max_idx = len(nuclei_index)
            valid_labels = props_df['label'] <= max_idx
            
            # Create mapping dict
            # label i corresponds to nuclei_index[i-1]
            mapper = {i: nuclei_index[i-1] for i in props_df.loc[valid_labels, 'label'] if i > 0}
            
            props_df['cell_id'] = props_df['label'].map(mapper)
        except Exception as e:
            print(f"      [WARN] Could not map implicit labels to IDs: {e}")
            return None

    # Final cleanup
    if 'cell_id' in props_df.columns:
        props_df = props_df.set_index('cell_id')
        props_df = props_df.drop(columns='label', errors='ignore')
        return props_df
    else:
        print("      [WARN] Failed to map labels to cell IDs.")
        return None


def main():
    parser = argparse.ArgumentParser(description="Filter nuclei, add polygons, and compute morphology.")
    
    # Defaults
    default_json = "/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json"
    default_geojson = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
    default_sdata = "/mnt/europa/valerio/data/zarr_store/samples"
    
    parser.add_argument("--json", default=default_json, help="Path to sample JSON")
    parser.add_argument("--geojson_dir", default=default_geojson, help="Path to GeoJSON dir")
    parser.add_argument("--sdata_dir", default=default_sdata, help="Path to sdata dir")
    parser.add_argument("--table_name", default='nuclei_counts_nop', help="Input table name")
    parser.add_argument("--table_name_new", default='final_table', help="Output table name")
    parser.add_argument("--nuclei_key", default='nuclei_boundaries', help="Nuclei key")
    
    args = parser.parse_args()

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

            # 2. Check Table
            if args.table_name not in sdata:
                print(f"  [SKIP] Table '{args.table_name}' not found.")
                continue
            
            table = sdata[args.table_name].copy()
            
            # 3. Get Coordinate System
            if not list(sdata.coordinate_systems):
                print("  [SKIP] No coord system.")
                continue
            coord_sys = list(sdata.coordinate_systems)[0]

            # 4. Initialize Columns
            for col in ['in_treatment', 'to_discard', 'in_tissue']:
                table.obs[col] = False

            # 5. Load GeoJSON
            geojson_file = os.path.join(args.geojson_dir, f"{sample_key}.geojson")
            polys_loaded = False
            
            if os.path.exists(geojson_file):
                try:
                    sample_poly = gpd.read_file(geojson_file)
                    if not sample_poly.empty:
                        sample_poly = sample_poly.set_crs(None, allow_override=True)
                        sample_poly_parsed = ShapesModel.parse(sample_poly, transformations={coord_sys: Identity()})
                        sdata["polys"] = sample_poly_parsed
                        polys_loaded = True
                    else:
                        print(f"  [WARN] GeoJSON empty: {sample_key}")
                except Exception as e:
                    print(f"  [ERROR] GeoJSON error: {e}")
            else:
                print(f"  [INFO] No GeoJSON found for {sample_key}")

            # 6. Spatial Join (if polys exist)
            if polys_loaded:
                try:
                    nuclei_filtering = sopa.spatial.sjoin(
                        sdata, args.nuclei_key, 'polys', how="inner", predicate="intersects", target_coordinate_system=coord_sys
                    )
                    
                    if not nuclei_filtering.empty:
                        nuclei_per_poly = nuclei_filtering.groupby('name').apply(lambda x: x.index.tolist())
                        
                        for source_key in nuclei_per_poly.index:
                            col_name = None
                            if "fibre_trattate" in source_key: col_name = "in_treatment"
                            elif "infiammazione" in source_key: col_name = "to_discard"
                            elif condition_id in source_key: col_name = "in_tissue"
                            
                            if col_name:
                                nuclei_list = nuclei_per_poly[source_key]
                                valid_indices = table.obs.index.intersection(nuclei_list)
                                if len(valid_indices) > 0:
                                    table.obs.loc[valid_indices, col_name] = True
                                    print(f"    Mapped '{source_key}' -> '{col_name}' ({len(valid_indices)})")
                except Exception as e:
                    print(f"  [ERROR] Sjoin failed: {e}")

            # 7. --- MORPHOLOGY STEP ---
            features_df = compute_morphology(sdata, args.nuclei_key, coord_sys)
            if features_df is not None:
                cols_to_add = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
                # Filter cols that actually exist in the result
                cols_to_add = [c for c in cols_to_add if c in features_df.columns]
                
                if cols_to_add:
                    # Join features to table.obs
                    # We drop existing columns if they already exist to avoid suffix issues
                    table.obs = table.obs.drop(columns=[c for c in cols_to_add if c in table.obs.columns], errors='ignore')
                    table.obs = table.obs.join(features_df[cols_to_add], how='left')
                    print(f"    Added morphology features: {cols_to_add}")

            # 8. Save
            sdata[args.table_name_new] = table
            sdata.write_element(args.table_name_new, overwrite=True)
            print(f"  [SUCCESS] Saved '{args.table_name_new}'")

if __name__ == "__main__":
    main()
