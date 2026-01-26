import pandas as pd
import os
import spatialdata as sd
from skimage import io
import numpy as np
from pathlib import Path
from py_scripts.utils.utils_fun import read_from_json

# info about the samples
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')
# directory of the sdata
spatialdata_dir = "/mnt/europa/valerio/data/zarr_store/samples"
# directory to save the images
images_dir = "/mnt/europa/valerio/HE_images/color_corrected/samples"
image_key = "full_image"

# --- Processing Loop ---
for blocco, samples in samples_dict.items():
    for sample, _ in samples.items():
        sample_id = f"{blocco}_{sample}.zarr"
        spatialdata_path = os.path.join(spatialdata_dir, sample_id)
        print(f"Processing {sample_id}...")
        try:
            # 1. Load the SpatialData object
            # (We only need the images, so we can try to lazily load if supported, 
            # but standard read_zarr is safest)
            sdata = sd.read_zarr(spatialdata_path)
            # 2. Check if image exists
            if image_key not in sdata.images:
                print(f"  -> Warning: Key '{image_key}' not found in {sample_id}. Skipping.")
                continue
            # 3. Extract the image data
            # sdata.images[key] returns a DataArray (often dask-backed)
            # Dims are typically (c, y, x) for SpatialData
            img = sdata[image_key]['scale0']['image'].data.compute()
            if img.shape[0] == 3: 
                img = np.moveaxis(img, 0, -1) # Convert (3, Y, X) -> (Y, X, 3)
            # 5. Define Output Filename
            out_name = f"{sample_id}.tif"
            out_path = os.path.join(images_dir, out_name)
            # 6. Save
            io.imsave(out_path, img, check_contrast=False, compression='zlib')
            print(f"  -> Saved to {out_path}")
            # Clean up memory
            del img
            del sdata
        except Exception as e:
            print(f"  -> Error processing {sample_id}: {e}")



print("Done.")


#single sample first try, it work

# first single sample try
# sdata = sd.read_zarr(f"{spatialdata_dir}/blocco1_sham.zarr")
# img = sdata[image_key]['scale0']['image'].data.compute()
# # 4. Convert to Numpy and Correct Dimensions
# # TiffFile usually expects (Y, X, C) for RGB images, or (C, Y, X).
# 
# # Optional: Transpose if you want YXC (Standard for many viewers like ImageJ/Napari defaults often vary, 
# # but YXC is standard for 'visual' RGB tiffs)
# # Check shape: if index 0 is 3 (RGB), move it to the end.
# if img.shape[0] == 3: 
#     img = np.moveaxis(img, 0, -1) # Convert (3, Y, X) -> (Y, X, 3)
# 
# # 5. Define Output Filename
# out_name = "blocco1_sham.tif"
# out_path = os.path.join(images_dir, out_name)
# 
# # 6. Save
# io.imsave(out_path, img, check_contrast=False, compression='zlib')
