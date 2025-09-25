#!/bin/bash

# -----------------------------------------------------------------------------
# This script runs the auto_segment_script.py on a predefined list of 
# .zarr directories to perform segmentation.
#
# To add more files to process, simply add their full paths to the 
# zarr_list array below.
# -----------------------------------------------------------------------------

echo "--- Starting segmentation process ---"
echo ""

# Get the directory where this bash script is located.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# --- FIX: Define the project's root directory ---
# This navigates three levels up from the script's directory.
PROJECT_ROOT=$( realpath "$SCRIPT_DIR/../../../" )

# Define the list of .zarr files to be processed.
zarr_list=(
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_c26foxO.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_c26STAT3.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26murf1.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26SMAD23.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_c26murf1.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_c26STAT3.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26foxO.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26SMAD23.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco4_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26murf1.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26SMAD23.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26STAT3.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26foxO.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26STAT3.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_sham.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26foxO.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26murf1.zarr"
"/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26SMAD23.zarr"
)

# Loop through each file path in the zarr_list array.
for zarr_file in "${zarr_list[@]}"; do
  echo "--> Processing target: $zarr_file"
  
  # --- FIX: Set PYTHONPATH before running the Python script ---
  # This tells Python to add the project root to its search path.
  PYTHONPATH=$PROJECT_ROOT python "$SCRIPT_DIR/samples_segmentation.py" "$zarr_file"
  
  echo "--> Finished target: $zarr_file"
  echo ""
done

echo "--- All files have been processed. ---"

# paths:
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco1_c26foxO.zarr"
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
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco5_c26STAT3.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26foxO.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26STAT3.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_c26.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco7_sham.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26foxO.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26murf1.zarr"
# "/mnt/europa/valerio/data/zarr_store/blocchi/blocco9_c26SMAD23.zarr"