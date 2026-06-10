#!/bin/bash

# Define the blocks you want to process
BLOCCHI=("blocco1" "blocco2" "blocco3" "blocco4" "blocco5" "blocco6" "blocco7" "blocco9")

# (Optional) Or dynamically extract them from your JSON if you have jq installed:
# BLOCCHI=$(jq -r 'keys[]' /mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json)

for BLOCCO in "${BLOCCHI[@]}"; do
    echo "Starting Python script for $BLOCCO..."
    
    # adata export 
    PYTHONPATH="$PWD" python py_scripts/R_export/sdata_to_adata_export.py "$BLOCCO"
    
    # Check if python exited with an error
    if [ $? -ne 0 ]; then
        echo "Warning: Python script failed for $BLOCCO. Continuing to next..."
    fi
done

echo "All blocks finished!"