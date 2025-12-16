import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import json
import os
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import py_scripts.utils.plots_fun as pf
import importlib
# After making changes to your_module.py, reload it
#importlib.reload(pf)

zarr_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples"
#zarr_path = "/mnt/europa/valerio/data/zarr_store/samples"
sdata = sd.read_zarr(f"{zarr_path}/blocco2_c26")

OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/spaceranger'
#OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/my_segmentation'


# categories = ['CD4 T cells', 'CD14+ Monocytes', 'B cells']
# colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
# 
# # Store in AnnData object
# adata.uns['louvain_colors'] = colors
colors = ["#0072B2", "#009E73", "#D55E00", "#E69F00", 
          "#CC79A7", "#56B4E9", "#F0E442", "#9B59B6", 
          "#0099B4", "#DDAA33", "#CC3311", "#44AA99", 
          "#AA4499", "#332288", "#A3E635"
]
# try to add custom colors to .uns
sdata['segmentation_counts'].uns['labels_singler_colors'] = colors

which_cells = ['MuSC', 'Immune_Cells', 'Myonuclei_NMJ', 'Myonuclei_IIb']
summary_df = pf.sample_cells_by_type(sdata, table_key = "segmentation_counts", which_cells = which_cells, n_samples = 5)

# solo bordi dei nuclei
pf.plot_single_nuclei(
  summary_df,
  sdata, 
  output_dir = f"{OUTPUT_DIR}/b2_c26_buff=300/borders", 
  table_key = "segmentation_counts",
  shapes_key = "blocco2_c26_nuclei_boundaries",
  image_key = "blocco2_c26_full_image",
  coord_sys = "blocco2",
  buffer = 300,
  figsize = (15,15),
  fill_alpha = 0,
  color_key = "None"
)

# nuclei colorati per tipo cellulare
pf.plot_single_nuclei(
  summary_df,
  sdata, 
  output_dir = f"{OUTPUT_DIR}/b2_c26_buff=300/cell_colors", 
  table_key = "segmentation_counts",
  shapes_key = "blocco2_c26_nuclei_boundaries",
  image_key = "blocco2_c26_full_image",
  coord_sys = "blocco2",
  buffer = 300,
  figsize = (15,15),
  fill_alpha = 1,
  color_key = "labels_singler"
)
