import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import geopandas as gpd
import pandas as pd
import skimage.io
import numpy as np
import json
import os
from spatialdata.models import ShapesModel, Labels2DModel
from spatialdata.transformations import Identity, get_transformation
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf

# ==========================================
# CONFIGURATION
# ==========================================

# Path to your JSON file containing the sample list
JSON_PATH = "/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json"  # <--- UPDATE THIS

# Directories
TIFF_DIR = "/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks"
ZARR_DIR = "/mnt/europa/valerio/data/zarr_store/samples"
OUTPUT_ADATA_DIR = "/mnt/europa/valerio/data/zarr_store/adatas/arivis_segmentation_tables"
EXCLUDED_SAMPLES = ["blocco1_sham", "blocco2_c26"]

# Ensure output directories exist
os.makedirs(OUTPUT_ADATA_DIR, exist_ok=True)

# Disable auto-save for SOPA to prevent partial writes during processing
sopa.settings.auto_save_on_disk = False


def process_sample(sample_key, block_name):
    """
    Process a single sample: Load TIFF mask -> Convert to Poly -> Aggregate -> Extract Features -> Save.
    """
    print(f"\n{'='*60}")
    print(f"Processing: {sample_key} (Block/System: {block_name})")
    print(f"{'='*60}")

    # 1. Define File Paths
    # Assumption: TIFF files are named "{sample_key}_finalprediction.tiff"
    tiff_path = os.path.join(TIFF_DIR, f"{sample_key}_finalprediction.tiff")
    zarr_path = os.path.join(ZARR_DIR, f"{sample_key}.zarr")
    adata_out_path = os.path.join(OUTPUT_ADATA_DIR, f"{sample_key}.h5ad")

    # 2. Validation Checks
    if not os.path.exists(tiff_path):
        print(f"❌ SKIPPING: TIFF not found at {tiff_path}")
        return False
    if not os.path.exists(zarr_path):
        print(f"❌ SKIPPING: Zarr not found at {zarr_path}")
        return False

    try:
        # 3. Load Data
        print(f"-> Loading TIFF...")
        nuclei_arivis = skimage.io.imread(tiff_path).astype(int)
        
        print(f"-> Loading Zarr SpatialData...")
        sdata = sd.read_zarr(zarr_path)

        # 4. Prepare Transformation
        # We assume the coordinate system we want matches 'block_name' or the first available one
        # Based on your snippet, we use sdata.coordinate_systems[0] to get the transform from 'full_image'
        coord_sys = sdata.coordinate_systems[0] 
        print(f"-> Using coordinate system: {coord_sys}")
        
        if 'full_image' not in sdata.images:
             print(f"❌ 'full_image' missing in Zarr. Skipping.")
             return False

        transf = get_transformation(sdata["full_image"], coord_sys)

        # 5. Add Labels to SpatialData
        print(f"-> Adding 'label_nuclei_arivis' to sdata...")
        nuclei_arivis_parsed = Labels2DModel.parse(
            data=nuclei_arivis, 
            transformations={block_name: transf}, 
            dims=("y", "x")
        )  
        sdata["label_nuclei_arivis"] = nuclei_arivis_parsed 

        # 6. Convert to Polygons (Using your custom function)
        print("-> Converting mask to polygons (sf.precise_to_polygons)...")
        nuclei_shapes = sf.precise_to_polygons(sdata["label_nuclei_arivis"])
        
        # Parse polys (Reset transformation to Identity for the vector layer)
        nuclei_shapes_parsed = ShapesModel.parse(
            nuclei_shapes, 
            transformations={block_name: Identity()}
        )
        sdata["nuclei_arivis_poly"] = nuclei_shapes_parsed

        # 7. SOPA Aggregation
        print("-> Running SOPA aggregation...")
        # Setup attributes
        sopa.utils.set_sopa_attrs(sdata, 
                    cell_segmentation_key='full_image', 
                    tissue_segmentation_key='full_image', 
                    transcripts_key=None, 
                    boundaries_key='nuclei_arivis_poly', 
                    bins_table_key='filtered'
        )

        # Run aggregation
        sopa.aggregate(sdata, key_added='arivis_nuclei', bins_key="filtered",
            shapes_key="nuclei_arivis_poly", expand_radius_ratio=0, min_transcripts=10,
            min_intensity_ratio=0.15, no_overlap=True
        )

        # 8. Feature Extraction (Morphology)
        print("-> Rasterizing polygons for feature extraction...")
        # Get extent using the correct coordinate system (block_name)
        element_extent = sd.get_extent(sdata['nuclei_arivis_poly'], coordinate_system=block_name, exact=True)
        
        # Rasterize back to image for regionprops
        sdata['raster_arivis_nuclei'] = sd.rasterize(
            sdata['nuclei_arivis_poly'],
            axes=["x", "y"],
            min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
            max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
            target_coordinate_system=block_name,
            target_unit_to_pixels=1,
        )

        label_mask = sdata['raster_arivis_nuclei'].values.squeeze().astype(np.int32)
        
        print("-> Calculating regionprops...")
        properties_to_extract = [
            'label', 'area', 'eccentricity', 'solidity', 'extent',
            'major_axis_length', 'minor_axis_length'
        ]
        props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))

        # Map labels to cell IDs
        if 'label_index_to_category' in sdata['raster_arivis_nuclei'].attrs:
            label_to_id = sdata['raster_arivis_nuclei'].attrs['label_index_to_category']
            props_df['cell_id'] = props_df['label'].map(label_to_id)
        else:
            props_df['cell_id'] = props_df['label']

        props_df = props_df.set_index('cell_id')
        props_df = props_df.drop(columns='label', errors='ignore')

        # Join to AnnData
        cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
        
        # Ensure index types match before joining
        if not props_df.empty:
            sdata['arivis_nuclei'].obs.index = sdata['arivis_nuclei'].obs.index.astype(props_df.index.dtype)
            sdata['arivis_nuclei'].obs = sdata['arivis_nuclei'].obs.join(props_df[cols], how='left')

        # 9. Save Updates to Zarr
        print("-> Saving elements to Zarr store (overwriting)...")
        # We explicitly delete and write if overwrite isn't behaving, but overwrite=True usually works
        sdata.write_element('label_nuclei_arivis')
        sdata.write_element('raster_arivis_nuclei')
        sdata.write_element('nuclei_arivis_poly')
        sdata.write_element('arivis_nuclei')

        # 10. Export AnnData
        print(f"-> Exporting AnnData to {adata_out_path}")
        adata = sdata['arivis_nuclei'].copy()
        
        # Add spatial coordinates to obs
        if 'spatial' in adata.obsm:
            adata.obs['x_coord'] = adata.obsm['spatial'][:, 0]
            adata.obs['y_coord'] = adata.obsm['spatial'][:, 1]
            del adata.obsm # Clean up
        
        adata.write_h5ad(adata_out_path)

        print("✅ SUCCESS")
        return True

    except Exception as e:
        print(f"❌ ERROR processing {sample_key}: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

# ==========================================
# MAIN EXECUTION
# ==========================================
# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == "__main__":
    
    if not os.path.exists(JSON_PATH):
        print(f"Critical Error: JSON file not found at {JSON_PATH}")
        exit(1)

    with open(JSON_PATH, 'r') as f:
        data_map = json.load(f)

    # Iterate through blocks and samples
    for block_name, samples_dict in data_map.items():
        for internal_name, details in samples_dict.items():
            
            # Extract the sample key (e.g., "blocco2_c26")
            sample_key = details.get("sample_key")
            
            if not sample_key:
                print(f"⚠️ Warning: No 'sample_key' found for {block_name} -> {internal_name}")
                continue

            # CHECK EXCLUSION LIST
            if sample_key in EXCLUDED_SAMPLES:
                print(f"⏩ SKIPPING excluded sample: {sample_key}")
                continue

            # Run processing
            process_sample(sample_key, block_name)

    print("\nBatch processing completed.")
