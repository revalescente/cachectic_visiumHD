import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# import spatialdata as sd

import tacco as tc

# plot settings
highres = True
default_dpi = 100.0
if highres:
    matplotlib.rcParams['figure.dpi'] = 648.0
    hr_ext = '_hd'
else:
    matplotlib.rcParams['figure.dpi'] = default_dpi
    hr_ext = ''

axsize = np.array([4,3])*0.5

region_colors = tc.pl.get_default_colors([f'region_{i}' for i in range(4)], offset=17)
split_names = np.array([f'sub_{i}' for i in range(4)])
split_colors = tc.pl.get_default_colors(split_names, offset=12)

# annotation colors
# The dictionary with the corrected key
color_dictionary = {
    "Endothelial": "#0072B2",
    "FAPs": "#009E73",
    "Immune_Cells": "#D55E00",
    "MuSC": "#E69F00",
    "Myonuclei_IIb": "#CC79A7",
    "Myonuclei_IIx": "#56B4E9",
    "Myonuclei_IIx_IIa": "#F0E442",
    "Myonuclei_IIx_IIb": "#9B59B6",
    "Myonuclei_MTJ": "#0099B4",
    "Myonuclei_NMJ": "#DDAA33",
    "Myonuclei_Trim63": "#CC3311",  # <-- CORRECTED: Changed "Trim6s3" to "Trim63"
    "Nervous_System": "#44AA99",
    "Pericyte": "#AA4499",
    "Smooth_Muscular": "#332288",
    "Tenocyte": "#A3E635"
}

# Create the pandas Series from the corrected dictionary
color_series = pd.Series(color_dictionary)

# Read reference and test data

adata_ref = ad.read_h5ad("/mnt/europa/valerio/data/multiome_sce/multiome_anndata.h5ad")
adata_test = ad.read_h5ad("/mnt/europa/valerio/data/zarr_store/adatas/samples_adata/blocco1_c26STAT3.h5ad")
# adata_temp = ad.read_h5ad("/mnt/europa/valerio/data/zarr_store/adatas/samples_adata/blocco1_c26foxO.h5ad")
  
# read spatial coords
spat_coords = pd.read_csv("/mnt/europa/valerio/data/data_tables/spatial_coords.csv", header=None)

start_row = 24166 # start from 0 so it's the row 24167
num_rows_to_select = 31755
end_row = start_row + num_rows_to_select

# Use .iloc to select rows by their integer position
adata_coords = spat_coords.iloc[start_row:end_row]

# --- Step 3: Add the NumPy array to adata.obsm ---
# The standard key for storing spatial coordinates is 'spatial'.
adata_test.obsm['spatial'] = adata_coords[[0, 1]].to_numpy()

# plot to control ----

sc.pl.spatial(
    adata_test,
    color="area",  # e.g., "leiden" or a gene name
    spot_size=40
)
plt.title("Spatial Plot with Nuclei Centroids")

plt.savefig("/mnt/europa/valerio//repositories/cachetic_visiumHD/figures/test_coords.png", dpi=300, bbox_inches='tight')


# CELLama not working 
# STEM 2 years without update in the repo - no env info
# Annotatability 11 months ago update - no env info
# TRY TACCO

# Annotate the spatial data with compositions of cell types -------------------------------------
 
adat_ref.varm['subclusters'] = pd.DataFrame(
    adata_ref.X.T.toarray().astype(np.float32), 
    index=adata_ref.var.index, 
    columns=adata_ref.obs['subclusters']
)

tc.tl.annotate(adata = adata_test, reference = adata_ref, annotation_key = 'subclusters', result_key='ClusterName')

# plot of subsample
adata_sub = adata_test[tc.sum(adata_test.X,axis=1)>=50].copy() # restrict downstream analysis to "good" beads
fig = tc.pl.scatter(adata_sub, keys='ClusterName', position_key = 'spatial', colors = color_series, joint=True, point_size=1, noticks=True, axes_labels=['X','Y']);

plt.savefig("/mnt/europa/valerio//repositories/cachetic_visiumHD/figures/test_tacco.png", dpi=300, bbox_inches='tight')



# cercava le info nel .varm (??) non ho una matrice di cluster quindi bo. ha comunque recuperato i subcluster giusti

cluster2type = adata_ref.obs[['subclusters','clusters']].drop_duplicates().groupby('clusters')['subclusters'].agg(lambda x: list(x.to_numpy()))
type2long = reference.obs[['type','long']].drop_duplicates().groupby('long')['type'].agg(lambda x: list(x.to_numpy()))

tc.utils.merge_annotation(adata_test, 'ClusterName', cluster2type, 'clusters')



# Annotate the spatial data with compositions of individual cells -------------------------------

# 
adata_ref.obs['cell'] = adata_ref.obs.index

adat_ref.varm['cell'] = pd.DataFrame(adata_ref.X.T.toarray().astype(np.float32), index=adata_ref.var.index, columns=adata_ref.var.index)
tc.tl.annotate(puck,adata_ref,'cell',result_key='cell',multi_center=None,)
