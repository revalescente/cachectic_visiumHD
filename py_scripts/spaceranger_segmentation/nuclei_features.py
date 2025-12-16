import spatialdata as sd
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf

# TO CHECK BEFORE RUNNING AGAIN, nuclei are in hires scale, not fullres! (so they are smaller, transform to coord sys fullres before use!)

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_sham")


nuclei_key = "blocco1_sham_nuclei_boundaries"
raster_key = "raster_nuclei"
blocco_key = "downscale_to_hires"

# features_df = sf.features_extraction(sdata, nuclei_element_name=nuclei_key)
print("Successfully extracted features.\n")


element_extent = sd.get_extent(sdata[nuclei_key], coordinate_system=blocco_key, exact=True)
sdata[raster_key] = sd.rasterize(
    sdata[nuclei_key],
    axes=["x", "y"],
    min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
    max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
    target_coordinate_system=blocco_key,
    target_unit_to_pixels=1,
)
# 4. --- Prepare Label Mask ---
# Squeeze the array to (H, W) and ensure it's an integer type
label_mask = sdata[raster_key].values.squeeze().astype(np.int32)
# 5. --- Compute Shape Features (regionprops) ---
properties_to_extract = [
    'label', 'area', 'eccentricity', 'solidity', 'extent',
    'major_axis_length', 'minor_axis_length'
]
props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))
# 6. --- Map Labels to Original Cell IDs ---

label_to_id = sdata[raster_key].attrs['label_index_to_category']
props_df['cell_id'] = props_df['label'].map(label_to_id)
   
props_df.to_csv('/mnt/europa/valerio/data/data_tables/b1_sham_nuclei_morph_SpaceRanger.csv', index=False)
