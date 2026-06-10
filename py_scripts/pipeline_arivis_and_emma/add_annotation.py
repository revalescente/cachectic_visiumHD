import pandas as pd
import spatialdata as sd
from py_scripts.utils.utils_fun import read_from_json

samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

nuclei_ann = pd.read_parquet("/mnt/europa/valerio/data/data_tables/arivis_plus_bins_annotation/arivis_nuclei_annotation.parquet")
bin_ann    = pd.read_parquet("/mnt/europa/valerio/data/data_tables/arivis_plus_bins_annotation/bin_016um_annotation.parquet")

# Set cell_id as index
nuclei_ann = nuclei_ann.set_index('nuclei_names')
bin_ann    = bin_ann.set_index('bin_names')

# ── CONFIG ──────────────────────────────────────────────────────────────────
sdata_dir    = "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins"
nuclei_table = "arivis_nuclei_table"
bin_table    = "square_016um"
# ────────────────────────────────────────────────────────────────────────────

for blocco, samples in samples_dict.items():
    for sample, _ in samples.items():
        # Load spatialdata object
        sdata = sd.read_zarr(f"{sdata_dir}/{blocco}_{sample}")
        
        # 1. Format keys for filtering correctly
        b_num = blocco.replace("blocco", "")
        nuc_sample_key = f"{blocco}_{sample}"
        bin_sample_key = f"{sample}_b{b_num}"  # Maps to format like 'c26STAT3_b1'
        
        # 2. Filter annotations using the correct formats
        sample_nuc = nuclei_ann[nuclei_ann["sample"] == nuc_sample_key][["nuclei_types"]]
        sample_bin = bin_ann[bin_ann["sample"] == bin_sample_key][["bin_types"]]
        
        # 3. Get the tables
        nuclei_counts = sdata.tables[nuclei_table]
        bin_counts    = sdata.tables[bin_table]
        
        # 4. FIX: Drop existing columns to avoid ValueError: columns overlap
        if "nuclei_types" in nuclei_counts.obs.columns:
            nuclei_counts.obs = nuclei_counts.obs.drop(columns=["nuclei_types"])
        if "bin_types" in bin_counts.obs.columns:
            bin_counts.obs = bin_counts.obs.drop(columns=["bin_types"])
        
        # 5. Join annotations into .obs
        nuclei_counts.obs = nuclei_counts.obs.join(sample_nuc, how="left")
        bin_counts.obs    = bin_counts.obs.join(sample_bin, how="left")
        
        # 6. Fill NaN for unmatched cells
        nuclei_counts.obs["nuclei_types"] = nuclei_counts.obs["nuclei_types"].fillna("Unknown")
        bin_counts.obs["bin_types"]       = bin_counts.obs["bin_types"].fillna("Unknown")
        
        # 7. Update and save back to zarr
        sdata.tables[nuclei_table] = nuclei_counts
        sdata.tables[bin_table]    = bin_counts
        
        sdata.delete_element_from_disk(nuclei_table)
        sdata.delete_element_from_disk(bin_table)
        sdata.write_element(nuclei_table)
        sdata.write_element(bin_table)
        
        # Clean up memory
        del sdata
        del sample_nuc
        del sample_bin
        
        print(f"Done! {blocco}_{sample} updated.")
        print("--- Nuclei Types ---")
        print(nuclei_counts.obs["nuclei_types"].value_counts())
        print("--- Bin Types ---")
        print(bin_counts.obs["bin_types"].value_counts())
        print("\n")
