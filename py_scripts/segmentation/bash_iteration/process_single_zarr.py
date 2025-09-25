#!/usr/bin/env python
# process_single_zarr.py

import os
import sys
import importlib
import spatialdata as sd
from sopa.io.standardize import read_zarr_standardized
import py_scripts.segmentation.segm_functions as sf

# --- A dictionary to map function names to actual, callable functions ---
# This is the key to making your script flexible.
# You can add any function you want to run on your sdata object here.
AVAILABLE_FUNCTIONS = {}

try:
    # Now that the path is set, this import should succeed.
    segm_module = importlib.import_module("py_scripts.segmentation.segm_functions")
    
    # Load the functions into the dictionary
    AVAILABLE_FUNCTIONS["postprocess_step"] = segm_module.postprocess_step
    AVAILABLE_FUNCTIONS["segmentation_step"] = segm_module.segmentation_step
    print(f"Successfully loaded functions: {list(AVAILABLE_FUNCTIONS.keys())}")

except ImportError as e:
    print(f"FATAL: Could not import function module: {e}")
    print("Please check that the file 'py_scripts/segmentation/segm_functions.py' exists relative to the script.")
    sys.exit(1)
except AttributeError as e:
    print(f"FATAL: A function was not found in the imported module: {e}")
    print("Please check the function names inside 'segm_functions.py'.")
    sys.exit(1)


def main():
    """
    Main function to read a single Zarr store, apply a named function, and save the result.
    """
    # This script now expects 2 arguments: a combined sample key and a function name.
    if len(sys.argv) != 3:
        print("Usage: python process_single_zarr.py <blocco_samplekey> <function_name>")
        print("Example: python process_single_zarr.py blocco1_c26 postprocess_step")
        sys.exit(1)

    combined_key = sys.argv[1]
    function_name = sys.argv[2]

    # --- Validate that the requested function is available ---
    if function_name not in AVAILABLE_FUNCTIONS:
        print(f"Error: Function '{function_name}' is not defined in AVAILABLE_FUNCTIONS.")
        print(f"Available options are: {list(AVAILABLE_FUNCTIONS.keys())}")
        sys.exit(1)
        
    function_to_run = AVAILABLE_FUNCTIONS[function_name]

    # --- Extract blocco and sample_key from the combined key ---
    try:
        # Split the string at the last underscore to handle cases like 'blocco_A_1_sample_B'
        blocco, sample_key = combined_key.rsplit('_', 1)
    except ValueError:
        print(f"Error: The combined key '{combined_key}' is not in the expected 'blocco_samplekey' format.")
        sys.exit(1)
        
    print(f"--- Starting Python session for: {combined_key} (PID: {os.getpid()}) ---")
    print(f"   - Blocco: {blocco}")
    print(f"   - Sample Key: {sample_key}")
    print(f"   - Applying function: {function_name}")

    # 1. Construct the path for this specific sample
    path_sdata = f"/mnt/europa/valerio/data/zarr_store/blocchi/{combined_key}.zarr"

    if not os.path.exists(path_sdata):
        print(f"Error: Data path not found at '{path_sdata}'")
        sys.exit(1)

    # 2. Read the data
    print(f"Reading data from: {path_sdata}")
    sdata = read_zarr_standardized(path_sdata)
    
    # 3. Apply a function of your choice to the spatial data, 
    # If it's from sopa it will be automatically saved
    # The result of this function is what gets saved.
    print(f"Executing '{function_name}'...")
    function_to_run(sdata) # This now correctly calls the function
    
    print(f"Successfully processed {blocco}_{sample_key}.")
    print(f"--- Python session for {blocco}_{sample_key} is closing. --- \n")


if __name__ == "__main__":
    main()
