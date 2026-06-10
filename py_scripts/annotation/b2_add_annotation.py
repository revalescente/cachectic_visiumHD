sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco2_c26")
# 16um annotation vector
annotations = pd.read_parquet("/mnt/europa/valerio/data/data_tables/arivis_plus_bins_annotation/bin_016um_annotation.parquet")

table_bin = sdata["square_016um"]
# 
annotations = annotations.set_index('bin_names')
annotations_c26_b2 = annotations[annotations['sample'] == 'c26_b2']
# Merge: add all annotation columns to . obs
table_bin.obs = table_bin.obs.join(annotations_c26_b2, how='left')
# Replace NaN with "Unknown" (or "NA", "") in the problematic column
table_bin.obs['bin_types'] = table_bin.obs['bin_types'].fillna('Unknown')
# Update and save
sdata.tables['square_016um'] = table_bin

sdata.delete_element_from_disk('square_016um')
sdata.write_element('square_016um')

# nuclei - manca il blocco 2 dei nuclei
annotations = pd.read_parquet("/mnt/europa/valerio/data/data_tables/arivis_plus_bins_annotation/arivis_nuclei_annotation.parquet")

table_nuc = sdata["arivis_nuclei_table"]
# 
annotations = annotations.set_index('nuclei_names')
annotations_c26_b2 = annotations[annotations['sample'] == 'blocco2_c26']
# Merge: add all annotation columns to . obs
table_nuc.obs = table_nuc.obs.join(annotations_c26_b2, how='left')
# Replace NaN with "Unknown" (or "NA", "") in the problematic column
table_nuc.obs['nuclei_types'] = table_nuc.obs['bin_types'].fillna("Unknown")
# Update and save
sdata.tables['arivis_nuclei_table'] = table_nuc

sdata.delete_element_from_disk('arivis_nuclei_table')
sdata.write_element('arivis_nuclei_table')
