import anndata
import spatialdata as sd
import numpy as np
import geopandas as gpd
import pandas as pd
import sopa
import json
import scanpy as sc 

sdata = sd.read_zarr('/mnt/europa/valerio/data/zarr_store/samples/blocco2_c26.zarr')
adata = sdata['nuclei_counts'].copy()

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-", "Mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL", "Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", "^Hb[^(p)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_genes(adata, min_cells=1)

# filter by nuclei area dimension
mask = sdata['nuclei_counts'].obs['area'] < 500
sdata['nuclei_counts'] = sdata['nuclei_counts'][mask]

#if not annotated
sdata.update_annotated_regions_metadata('nuclei_counts', region_key='region')
sdata.set_table_annotates_spatialelement('nuclei_counts', region = 'filtered_nuclei', region_key='region', instance_key='cell_id')

# saving new objects
nuclei_boundaries = sdata['nuclei_boundaries'].copy()
nuclei_boundaries.reset_index(inplace=True)
nuclei_boundaries.rename(columns={'index': 'cell_id'}, inplace=True)
cell_ids_to_keep = set(sdata['nuclei_counts'].obs['cell_id'])
filtered_boundaries = nuclei_boundaries[nuclei_boundaries['cell_id'].isin(cell_ids_to_keep)]

tmp_table = sdata['nuclei_counts'].copy()

sdata.delete_element_from_disk('filtered_nuclei')
sdata.delete_element_from_disk('nuclei_counts')

sdata['filtered_nuclei'] = filtered_boundaries
sdata['nuclei_counts'] = tmp_table

sdata.write_element('filtered_nuclei', overwrite=False, format=None)
sdata.write_element('nuclei_counts', overwrite=False, format=None)
