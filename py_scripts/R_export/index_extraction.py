import spatialdata as sd
import numpy as np
import pandas as pd
from py_scripts.utils.utils_fun import read_from_json

samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
# data
# sdata = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys}.zarr")
results_list = []

for block, samples in samples_dict.items():
    for condition, details in samples.items():
        sample_key = details['sample_key']
        zarr_path = f"/mnt/europa/valerio/data/zarr_store/samples/{sample_key}.zarr"
        
        print(f"Processing: {sample_key} ...")
        
        try:
            sdata = sd.read_zarr(zarr_path)
            
            if 'nuclei_counts' in sdata.tables:
                adata = sdata.tables['nuclei_counts']
                obs = adata.obs
                
                # 2. Create a DataFrame for this sample
                # We initialize it with the Cell IDs and the Sample Key
                df_sample = pd.DataFrame({
                    'sample_key': sample_key,
                    'cell_id': obs.index.tolist()
                })
                
                # 3. Extract 'in_treatment' (Assuming it exists, or filling False if missing)
                if 'in_treatment' in obs.columns:
                    df_sample['in_treatment'] = obs['in_treatment'].values
                else:
                    df_sample['in_treatment'] = False # Fallback
                
                # 4. Extract 'GFP_value' with conditional logic
                if 'GFP_value' in obs.columns:
                    df_sample['GFP_value'] = obs['GFP_value'].values
                else:
                    # Create a column of False with the same length as the dataframe
                    df_sample['GFP_value'] = False
                
                # Add to our list
                results_list.append(df_sample)
                print(f"  -> Extracted {len(df_sample)} cells.")
                
            else:
                print(f"  -> WARNING: 'nuclei_counts' table not found in {sample_key}")
            
            del sdata
            
        except Exception as e:
            print(f"  -> ERROR reading {sample_key}: {e}")

# 5. Concatenate everything into one big table
final_df = pd.concat(results_list, ignore_index=True)

# 6. Save to CSV for R
output_csv = "/mnt/europa/valerio/data/data_tables/nuclei_indices_after_GFP_poly_filter.csv"
final_df.to_csv(output_csv, index=False)
print(f"Saved full metadata to {output_csv}")


# save the indices

df_export = pd.DataFrame(
    [(sample, cid) for sample, ids in cell_indices.items() for cid in ids],
    columns=['sample_key', 'cell_id']
)

# Save to CSV
df_export.to_csv("/mnt/europa/valerio/data/data_tables/nuclei_indices_after_GFP_poly_filter.csv", index=False)
print(f"Saved {len(df_export)} rows to CSV.")
