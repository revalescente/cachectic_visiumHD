import pandas as pd
import os
import spatialdata as sd
from py_scripts.utils.utils_fun import read_from_json


samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

#spatialdata_dir = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples"
spatialdata_dir = "/mnt/europa/valerio/data/zarr_store/samples"
#table_key = "segmentation_counts"  # Replace with the table key of your choice
table_key = "nuclei_counts_nop"  # Replace with the table key of your choice
new_table_key = "nuclei_counts"
#parquet_dir = "/mnt/europa/valerio/data/data_tables/spaceranger/nuclei_filtering_info"
parquet_dir = "/mnt/europa/valerio/data/data_tables/my_segmentation/nuclei_filtering_info"

for blocco, samples in samples_dict.items():
  for sample,_ in samples.items():
    # sdata
    spatialdata_path = f"{spatialdata_dir}/{blocco}_{sample}.zarr"
    sdata = sd.read_zarr(spatialdata_path)
    # parquet
    parquet_path = f"{parquet_dir}/{blocco}_{sample}.parquet"
    new_cols_df = pd.read_parquet(parquet_path)
    if "cell_id" in new_cols_df.columns:
        new_cols_df = new_cols_df.set_index("cell_id")
    cols_to_use = new_cols_df.columns.difference(sdata[table_key].obs.columns)
    # 1. Join
    sdata[table_key].obs = sdata[table_key].obs.join(new_cols_df[cols_to_use], how='left')
    # 2. Strict Filter (Remove cells based on primary label)
    if 'labels_singler' in sdata[table_key].obs.columns:
        sdata[new_table_key] = sdata[table_key][sdata[table_key].obs['labels_singler'].notna()].copy()
    # 3. Safe Fill (Fix remaining string columns to prevent Zarr crash)
    # iterate over ONLY the new columns to avoid messing up existing data
    for col in cols_to_use:
        # Ensure column still exists (it might have been dropped during filtering)
        if col in sdata[new_table_key].obs.columns:
            # Check the dtype. We generally want to fix 'object', 'bool', and 'category'.
            # We skip 'int' or 'float' to avoid turning measurements into strings.
            if sdata[new_table_key].obs[col].dtype.kind in ['O', 'b', 'S', 'U']: # O=Object, b=Bool
                # 1. Fill NaNs with "NA" (or "False" if you prefer)
                # 2. Force the whole column to string. 
                #    This converts True -> "True", False -> "False", NaN -> "NA"
                sdata[new_table_key].obs[col] = sdata[new_table_key].obs[col].fillna("NA").astype(str)
    sdata.delete_element_from_disk(new_table_key)
    sdata.write_element(new_table_key)
    del sdata
    del new_cols_df


sdata = sd.read_zarr(f"{spatialdata_dir}/blocco7_c26.zarr")
df_try = pd.read_parquet(f"{parquet_dir}/blocco2_c26.parquet")
# parquet
if "cell_id" in df_try.columns:
    df_try = df_try.set_index("cell_id")
cols_to_use = df_try.columns.difference(sdata[table_key].obs.columns)
# 1. Join
sdata[table_key].obs = sdata[table_key].obs.join(df_try[cols_to_use], how='left')
# 2. Strict Filter (Remove cells based on primary label)
if 'labels_singler' in sdata[table_key].obs.columns:
    sdata[new_table_key] = sdata[table_key][sdata[table_key].obs['labels_singler'].notna()].copy()
# 3. Safe Fill (Fix remaining string columns to prevent Zarr crash)
# iterate over ONLY the new columns to avoid messing up existing data
for col in cols_to_use:
    # Ensure column still exists (it might have been dropped during filtering)
    if col in sdata[new_table_key].obs.columns:
        # Check the dtype. We generally want to fix 'object', 'bool', and 'category'.
        # We skip 'int' or 'float' to avoid turning measurements into strings.
        if sdata[new_table_key].obs[col].dtype.kind in ['O', 'b', 'S', 'U']: # O=Object, b=Bool
            # 1. Fill NaNs with "NA" (or "False" if you prefer)
            # 2. Force the whole column to string. 
            #    This converts True -> "True", False -> "False", NaN -> "NA"
            sdata[new_table_key].obs[col] = sdata[new_table_key].obs[col].fillna("NA").astype(str)

sdata.delete_element_from_disk(new_table_key)
sdata.write_element(new_table_key)
