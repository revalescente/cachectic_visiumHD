#!/bin/bash

# Define the blocks you want to process
BLOCCHI=("blocco1" "blocco3" "blocco4" "blocco5" "blocco6" "blocco7" "blocco9")

# (Optional) Or dynamically extract them from your JSON if you have jq installed:
# BLOCCHI=$(jq -r 'keys[]' /mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json)

for BLOCCO in "${BLOCCHI[@]}"; do
    echo "Starting Python script for $BLOCCO..."
    
    # preprocessing all samples with bins and arivis nulcei
    # PYTHONPATH="$PWD" python py_scripts/pipeline_arivis_and_emma/process_single_blocco.py "$BLOCCO"
    
    # process only the nuclei of arivis to add the spatial info (to discard and GFP region and values)
    PYTHONPATH="$PWD" python py_scripts/pipeline_arivis_and_emma/GFP_info_arivis_nuclei.py "$BLOCCO"
    
    # Check if python exited with an error
    if [ $? -ne 0 ]; then
        echo "Warning: Python script failed for $BLOCCO. Continuing to next..."
    fi
done

echo "All blocks finished!"

