import os
import spatialdata as sd
from skimage import io
import numpy as np

# Define the directories
input_dir = "/mnt/europa/valerio/data/zarr_store/binned_samples/version_1.0.0/"
output_dir = "/mnt/europa/valerio/data/zarr_store/binned_samples/version_1.0.1/"
tif_dir = "/mnt/europa/valerio/Fluo_images/warped_tif/samples"

# Ensure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(tif_dir, exist_ok=True)

# Iterate through the items in the input directory
for item in os.listdir(input_dir):
    input_zarr_path = os.path.join(input_dir, item)
    output_zarr_path = os.path.join(output_dir, item)
    
    # Check if it's a directory (Zarr stores are directories)
    if not os.path.isdir(input_zarr_path):
        continue
        
    print(f"\nProcessing: {item}")
    
    # --- NEW CHECK: Skip if already exists in output_dir ---
    if os.path.exists(output_zarr_path):
        print(f"  Skipping '{item}': already exists in {output_dir}")
        continue
    
    # 1. Read the SpatialData Zarr store lazily
    try:
        sdata = sd.read_zarr(input_zarr_path)
    except Exception as e:
        print(f"  [Error] Failed to read {item} as SpatialData Zarr: {e}")
        continue
        
    # 2. Save it back to the new version 1.0.1 folder
    print(f"  Saving Zarr to: {output_zarr_path} ...")
    try:
        # write() will save the exact same sdata structure to the new path
        sdata.write(output_zarr_path)
    except Exception as e:
        print(f"  [Error] Failed to write {item}: {e}")
        continue
        
    # 3. Check for fluo_image and save the highest resolution as .tif
    if 'fluo_image' in sdata.images:
        print(f"  'fluo_image' found. Extracting high-res and saving as .tif ...")
        try:
            # sdata.images[key] returns a DataArray (often dask-backed)
            # Dims are typically (c, y, x) for SpatialData
            img = sdata['fluo_image']['scale0']['image'].data.compute()
            if img.shape[0] == 3: 
                img = np.moveaxis(img, 0, -1) # Convert (3, Y, X) -> (Y, X, 3)

            # Save the numpy array as a TIF file
            tif_path = os.path.join(tif_dir, f"{item}_fluo_image.tif")
            io.imsave(tif_path, img, check_contrast=False, compression='zlib')
            print(f"  -> Saved to {output_zarr_path}")
            print(f"  Saved TIF to: {tif_path}")
            # Clean up memory
            del img
            del sdata

        except Exception as e:
            print(f"  [Error] Failed to process/save fluo_image for {item}: {e}")
    else:
        print(f"  No 'fluo_image' found in {item}, skipping TIF generation.")
        
print("\nProcessing complete!")