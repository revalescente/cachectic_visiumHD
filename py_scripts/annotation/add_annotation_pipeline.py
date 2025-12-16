import pandas as pd
import spatialdata as sd
import json

# Spaceranger qui --------------------------------------------------------------
#
#sotto my segmentation

# dictionary to manage the various samples
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    samples_dict = json.load(f)

# just blocco1 and blocco2
subset_dict = {k: samples_dict[k] for k in ['blocco1', 'blocco2'] if k in samples_dict}

for blocco, samples in subset_dict.items():
  for sample,_ in samples.items():
    sdata = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{blocco}_{sample}")
    annotations = pd.read_parquet(f"/mnt/europa/valerio/data/data_tables/annotation_results_filtered/spaceranger/{blocco}_{sample}.parquet")
    # Set cell_id as index
    annotations = annotations.set_index('cell_id')
    # Get nuclei_counts and set cell_id as index if needed
    nuclei_counts = sdata.tables['segmentation_counts']
    # Merge: add all annotation columns to . obs
    nuclei_counts.obs = nuclei_counts.obs.join(annotations, how='left')
    # Replace NaN with "Unknown" (or "NA", "") in the problematic column
    nuclei_counts.obs['labels_singler'] = nuclei_counts.obs['labels_singler'].fillna("Unknown")
    # Update and save
    sdata.tables['segmentation_counts'] = nuclei_counts
    sdata.delete_element_from_disk('segmentation_counts')
    sdata.write_element('segmentation_counts')
    del sdata
    del annotations


sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco1_c26foxO.zarr")
annotations = pd.read_parquet("/mnt/europa/valerio/data/data_tables/annotation_results_filtered/blocco1_c26foxO.parquet")
# Set cell_id as index
annotations = annotations.set_index('cell_id')

# Get nuclei_counts and set cell_id as index if needed
nuclei_counts = sdata.tables['nuclei_counts']

# Merge: add all annotation columns to . obs
nuclei_counts.obs = nuclei_counts.obs.join(annotations, how='left')
# Replace NaN with "Unknown" (or "NA", "") in the problematic column
nuclei_counts.obs['labels_singler'] = nuclei_counts.obs['labels_singler'].fillna("Unknown")

# Update and save
sdata.tables['nuclei_counts'] = nuclei_counts
sdata.delete_element_from_disk('nuclei_counts')
sdata.write_element('nuclei_counts')


# My segmentation ------------------------------------------------------------------

# dictionary to manage the various samples
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    samples_dict = json.load(f)
    
for blocco, samples in samples_dict.items():
  for sample,_ in samples.items():
    sdata = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{blocco}_{sample}.zarr")
    annotations = pd.read_parquet(f"/mnt/europa/valerio/data/data_tables/annotation_results_filtered/my_segmentation/{blocco}_{sample}.parquet")
    # Set cell_id as index
    annotations = annotations.set_index('cell_id')
    
    # Get nuclei_counts and set cell_id as index if needed
    nuclei_counts = sdata.tables['nuclei_counts']
    
    # Merge: add all annotation columns to . obs
    nuclei_counts.obs = nuclei_counts.obs.join(annotations, how='left')
    # Replace NaN with "Unknown" (or "NA", "") in the problematic column
    nuclei_counts.obs['labels_singler'] = nuclei_counts.obs['labels_singler'].fillna("Unknown")
    
    # Update and save
    sdata.tables['nuclei_counts'] = nuclei_counts
    sdata.delete_element_from_disk('nuclei_counts')
    sdata.write_element('nuclei_counts')
  
    del sdata
    del annotations


sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco1_c26foxO.zarr")
annotations = pd.read_parquet("/mnt/europa/valerio/data/data_tables/annotation_results_filtered/blocco1_c26foxO.parquet")
# Set cell_id as index
annotations = annotations.set_index('cell_id')

# Get nuclei_counts and set cell_id as index if needed
nuclei_counts = sdata.tables['nuclei_counts']

# Merge: add all annotation columns to . obs
nuclei_counts.obs = nuclei_counts.obs.join(annotations, how='left')
# Replace NaN with "Unknown" (or "NA", "") in the problematic column
nuclei_counts.obs['labels_singler'] = nuclei_counts.obs['labels_singler'].fillna("Unknown")

# Update and save
sdata.tables['nuclei_counts'] = nuclei_counts
sdata.delete_element_from_disk('nuclei_counts')
sdata.write_element('nuclei_counts')
