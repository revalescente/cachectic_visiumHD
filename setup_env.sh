#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
# The name of your environment YAML file.
ENV_YAML="deepl-segment-bio.yml"
# The name of the conda environment defined inside the YAML file.
# IMPORTANT: This MUST match the 'name:' field in your .yml file.
ENV_NAME="deepl-segment-bio"
# The path to your EXISTING local clone of the sopa repository.
# Change this path if your 'sopa' folder is not in the same directory as this script.
SOPA_DIR_PATH="/mnt/europa/valerio/repositories/Sopa"

# --- Main Script ---

echo "--- Step 1: Configuring Conda Channels ---"
# Prioritizing conda-forge is best practice for robust dependency resolution.
conda config --add channels conda-forge
conda config --set channel_priority strict
echo "Conda channels configured to prioritize 'conda-forge'."
echo

echo "--- Step 2: Creating Conda Environment '$ENV_NAME' from '$ENV_YAML' ---"
# Check if the environment already exists to avoid errors.
if conda env list | grep -q -w "$ENV_NAME"; then
    echo "Warning: Environment '$ENV_NAME' already exists. Skipping creation."
    echo "If you want to recreate it, first run: conda env remove --name $ENV_NAME"
else
    if [ ! -f "$ENV_YAML" ]; then
        echo "Error: The environment file '$ENV_YAML' was not found."
        echo "Please make sure the file is in the same directory as this script."
        exit 1
    fi
    conda env create -f "$ENV_YAML"
    echo "Success: Environment '$ENV_NAME' created successfully."
fi
echo

echo "--- Step 3: Installing sopa in Editable Mode ---"
# Check that the specified sopa directory actually exists before trying to install it.
if [ ! -d "$SOPA_DIR_PATH" ]; then
    echo "Error: The specified sopa directory '$SOPA_DIR_PATH' does not exist."
    echo "Please update the SOPA_DIR_PATH variable in the script to point to your cloned repository."
    exit 1
fi

# We use 'conda run' to execute the command inside the specific environment
# without needing to activate it in the script, which is more reliable.
conda run -n "$ENV_NAME" pip install -e "./$SOPA_DIR_PATH"
echo
echo "Success: Sopa has been installed in editable mode in the '$ENV_NAME' environment."
echo
echo "To activate your new environment, run:"
echo "conda activate $ENV_NAME"