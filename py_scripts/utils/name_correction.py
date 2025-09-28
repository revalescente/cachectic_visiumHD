import spatialdata as sd
import spatialdata_plot
from spatialdata import SpatialData
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import (Identity, get_transformation, set_transformation, remove_transformations_to_coordinate_system)
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import gc
import re
import sopa
import json
import time
import sys
from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf

# after reading the dataset i need to change the names of spatial elements - 
# then i concatenate and then i save so i don't need to save the intermediate step

# for the subsequent analysis I need:

# spatial elements
# - bloccokey_full_image

# - bloccokey_filtered_bins
# - bloccokey_filtered_nuclei
# - bloccokey_intissue
# - bloccokey_intissue_002um

# table
# - nuclei_counts_nop

# sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_c26foxO.zarr")

def unique_naming(sdata):
  
  # new suffix for the elements
  sample_key = sdata.path.stem
  
  # Use regex to extract blocco_key and samples_key
  match = re.match(r'(blocco\d+)_(\w+)', sample_key)
  if match:
      blocco, sampletype = match.group(1), match.group(2)
      
  else:
      raise ValueError(f"Could not parse blocco from: {sdata_path}")
      
  new_names = [
    "nuclei_counts_nop",
    f"{sample_key}_full_image",
    f"{sample_key}_filtered_bins",
    f"{sample_key}_filtered_nuclei",
    f"{sample_key}_intissue",
    f"{sample_key}_intissue_002um"
  ]
  
  old_names = [
    "nope",
    f"{blocco}_full_image",
    f"{blocco}_filtered_bins",
    f"{blocco}_filtered_nuclei",
    f"{blocco}_intissue",
    f"{blocco}_intissue_002um"
  ] 
    
  for old_name, new_name in zip(old_names[1:], new_names[1:]):
    # Check if the name contains '_full_image' to identify it as an image
    if "_full_image" in new_name:
      print(f"  Copying image: '{old_name}' -> '{new_name}'")
      sdata[new_name] = sdata[old_name].copy()
    else:
      print(f"  Deep copying element: '{old_name}' -> '{new_name}'")
      sdata[new_name] = sd.deepcopy(sdata[old_name])
      
  # Let's annotate the table of the genes vs nuclei with the filtered nuclei
  sdata["nuclei_counts_nop"].obs["region"] = f"{sample_key}_filtered_nuclei"
  sdata.set_table_annotates_spatialelement("nuclei_counts_nop", region = f"{sample_key}_filtered_nuclei", region_key="region", instance_key="cell_id")
  
  # matching table with the filtered nuclei
  sdata['nuclei_counts_nop'] = sd.match_table_to_element(sdata, element_name = f"{sample_key}_filtered_nuclei", table_name='nuclei_counts_nop')
    
  sdata_new = sdata.subset(new_names[1:], filter_tables = True, include_orphan_tables=False)
  
  # add some useful column: blocco, sample_key, and sample 
  sdata_new['nuclei_counts_nop'].obs['sample_id'] = sample_key
  sdata_new['nuclei_counts_nop'].obs['slide'] = f'{sample_key}_full_image'
  sdata_new['nuclei_counts_nop'].obs['blocco'] = f'{blocco}'
  sdata_new['nuclei_counts_nop'].obs['sampletype'] = f'{sampletype}'
  
  # add nuclei features to anndata.obs
  features_df = sf.features_extraction(sdata_new, nuclei_element_name = f"{sample_key}_filtered_nuclei")
  cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
  sdata_new['nuclei_counts_nop'].obs = sdata_new['nuclei_counts_nop'].obs.join(features_df[cols], how='left')
  
  # create the new coordinate_system called as the sample_key
  transformed = get_transformation(sdata_new.images[f'{sample_key}_full_image'], to_coordinate_system=f'{blocco}')
  for element in new_names[1:]:
    set_transformation(sdata_new[element], transformed, to_coordinate_system=f'{sample_key}')

  remove_transformations_to_coordinate_system(sdata_new, f"{blocco}")
  
  print(sdata_new)
  savepath = f'/mnt/europa/valerio/data/zarr_store/samples/{sample_key}.zarr'
  sdata_new.write(savepath)



def main():
    # Check if a command-line argument (the path) was provided
    if len(sys.argv) < 2:
        print("Error: Please provide the path to the .zarr file.")
        sys.exit(1)

    # The first argument is the script name, the second is our path
    zarr_path = sys.argv[1]
    
    print(f"Reading sdata from: {zarr_path}")
    
    try:
        # Read the .zarr file into a SpatialData object
        sdata = sd.read_zarr(zarr_path)
        
        # Call your processing function
        sdata = unique_naming(sdata)
        
        # Here you would typically save the modified sdata object, for example:
        # print(f"Saving modified sdata back to: {zarr_path}")
        # sdata.write(zarr_path, overwrite=True)
        
    except Exception as e:
        print(f"An error occurred while processing {zarr_path}: {e}")

if __name__ == "__main__":
    main()
