#!/bin/bash

# -----------------------------------------------------------------------------
# This script runs a Python module on a list of .zarr files.
# It ensures the Python script is run from the project root to resolve imports.
# -----------------------------------------------------------------------------

echo "--- Starting batch processing ---"
echo ""

# 1. Get the directory where this bash script is located.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# 2. Define the project's root directory by navigating up from the script's location.
#    (e.g., utils -> py_scripts -> project_root)
PROJECT_ROOT=$( realpath "$SCRIPT_DIR/../../" )

# Define the array of .zarr file paths to be processed.
path_list=(
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_c26foxO.zarr" # DONE:
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_c26STAT3.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_sham.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26murf1.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26SMAD23.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_c26murf1.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_c26STAT3.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_sham.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26foxO.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26SMAD23.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26murf1.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26SMAD23.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26STAT3.zarr" # DONE!
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26foxO.zarr" # PROBLEM
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26STAT3.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26foxO.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26murf1.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26SMAD23.zarr"
)

# 3. Change the current directory to the project root.
#    This is the crucial step that makes `python -m` work.
echo "Changing directory to project root: $PROJECT_ROOT"
cd "$PROJECT_ROOT"

# Check if the cd was successful
if [ $? -ne 0 ]; then
  echo "Error: Could not change to project root directory. Aborting."
  exit 1
fi

# 4. Loop through each path and call the Python script AS A MODULE.
for zarr_path in "${path_list[@]}"; do
  echo "--> Processing target: $zarr_path"
  
  # Use `python -m` with the module path (dots instead of slashes, no .py)
  # The script name_correction.py is in py_scripts/utils/
  python -m py_scripts.utils.name_correction "$zarr_path"
  
  echo "--> Finished target: $zarr_path"
  echo ""
done

echo "--- All files have been processed. ---"