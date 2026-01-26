import spatialdata as sd
import anndata as an
import py_scripts.utils.spaceranger_utility as su
import os
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.segmentation.segm_functions as sf


samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')


spatialdata_dir = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples"
table_key = "segmentation_counts"  # Replace with the table key of your choice
output_dir = "/mnt/europa/valerio/data/zarr_store/adatas/spaceranger_tables"

for blocco, samples in samples_dict.items():
  for sample,_ in samples.items():
    spatialdata_path = f"{spatialdata_dir}/{blocco}_{sample}"
    sdata = sd.read_zarr(spatialdata_path)
    prop_df = sf.features_extraction(sdata, nuclei_element_name = f"{blocco}_{sample}_nuclei_boundaries")
    # Identify which columns are missing in the .obs DataFrame
    missing_columns = [col for col in prop_df.columns if col not in sdata[table_key].obs.columns]
    # Add only the missing columns to .obs
    for col in missing_columns:
        if col in prop_df:  # Ensure the column also exists in df
            sdata[table_key].obs[col] = prop_df[col]
    sdata.delete_element_from_disk('segmentation_counts')
    sdata.write_element('segmentation_counts')
    del sdata
    su.save_table_as_h5ad(spatialdata_path, table_key, output_dir)



spatialdata_path = f"{spatialdata_dir}/blocco1_c26STAT3"
sdata = sd.read_zarr(spatialdata_path)
prop_df = sf.features_extraction(sdata, nuclei_element_name = "blocco1_c26STAT3_nuclei_boundaries")
# Identify which columns are missing in the .obs DataFrame
missing_columns = [col for col in prop_df.columns if col not in sdata[table_key].obs.columns]
# Add only the missing columns to .obs
for col in missing_columns:
    if col in prop_df:  # Ensure the column also exists in df
        sdata[table_key].obs[col] = prop_df[col]

sdata.delete_element_from_disk('segmentation_counts')
sdata.write_element('segmentation_counts')
su.save_table_as_h5ad(spatialdata_path, table_key, output_dir)
