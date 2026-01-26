import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import json
import os
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import py_scripts.utils.plots_fun as pf
import sopa
from importlib import reload
#reload(pf)

#zarr_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples"
zarr_path = "/mnt/europa/valerio/data/zarr_store/samples"
sdata = sd.read_zarr(f"{zarr_path}/blocco2_c26.zarr")

#OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/spaceranger'
OUTPUT_DIR = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/sampled_cells_plots/my_segmentation'

# COORD_SYS = 'blocco2'
# TABLE_KEY = 'segmentation_counts'
# IMAGE_KEY = f'{sdata.path.stem}_full_image'
# SHAPES_KEY = f'{sdata.path.stem}_nuclei_boundaries'

COORD_SYS = 'blocco2'           # From your sdata printout
TABLE_KEY = 'nuclei_counts'
IMAGE_KEY = 'full_image'        # The background image
SHAPES_KEY = 'nuclei_boundaries'  # The shapes matching your table

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
sdata[TABLE_KEY].uns['labels_singler_colors'] = colors

which_cells = ['MuSC', 'Immune_Cells', 'Myonuclei_NMJ', 'Myonuclei_IIb']
summary_df = pf.sample_cells_by_type(sdata, table_key = TABLE_KEY, which_cells = which_cells, n_samples = 5)

# solo bordi dei nuclei
pf.plot_single_nuclei(
  summary_df,
  sdata, 
  output_dir = f"{OUTPUT_DIR}/{sdata.path.stem}_buff=100/borders", 
  table_key = TABLE_KEY,
  shapes_key = SHAPES_KEY,
  image_key = IMAGE_KEY,
  coord_sys = COORD_SYS,
  buffer = 100,
  figsize = (15,15),
  fill_alpha = 0,
  color_key = "None"
)

# nuclei colorati per tipo cellulare
pf.plot_single_nuclei(
  summary_df,
  sdata, 
  output_dir = f"{OUTPUT_DIR}/{sdata.path.stem}_buff=100/cell_colors", 
  table_key = TABLE_KEY,
  shapes_key = SHAPES_KEY,
  image_key = IMAGE_KEY,
  coord_sys = COORD_SYS,
  buffer = 100,
  figsize = (15,15),
  fill_alpha = 1,
  color_key = "labels_singler"
)

# ---------

# filter out to_discard

my segmentation blocco2 da sistemare?
sdata["intissue_poly"] = sdata["GFP_poly"][sdata["GFP_poly"].name == "c26"]
intissue_nuclei = sopa.spatial.sjoin(
            sdata,
            "nuclei_boundaries",
            "intissue_poly",
            how='inner',
            predicate='intersects',
            target_coordinate_system=COORD_SYS
        )
sdata["nuclei_boundaries"] = intissue_nuclei.copy()
sdata["nuclei_counts"] = sd.match_table_to_element(sdata, "nuclei_boundaries", table_name='nuclei_counts')

keep_mask = sdata[TABLE_KEY].obs['to_discard'].astype(str) == "False"
print(keep_mask.value_counts())
sdata["filtered"] = sdata[TABLE_KEY][keep_mask].copy()
print(f"Remaining cells: {sdata['filtered'].n_obs}")

keep_mask = (sdata[TABLE_KEY].obs["area"] > 50) & (sdata[TABLE_KEY].obs["area"] < 1001)
sdata["filtered"] = sdata[TABLE_KEY][keep_mask].copy()


# sr no problem instead.... bo
keep_mask = sdata[TABLE_KEY].obs['to_discard'] == False
sdata["filtered"] = sdata[TABLE_KEY][keep_mask].copy()

# Various QC plots with black background to improve visibility

fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="AspectRatio",
    table_name="filtered",
    cmap="plasma"
).pl.show(
    ax=ax, 
    title=f"Aspect Ratio {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/AspectRatio_{sdata.path.stem}.png"
plt.savefig(save_path, facecolor='black', dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="subsets_Mito_percent",
    table_name="filtered",
    cmap="plasma"
).pl.show(
    ax=ax, 
    title=f"Mitocondrial genes percentage {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/Mito_percent_{sdata.path.stem}.png"
plt.savefig(save_path, facecolor='black', dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="log2CountArea",
    table_name="filtered",
    cmap="plasma"
).pl.show(
    ax=ax, 
    title=f"log2 ratio between Counts and Area {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/log2CountArea_{sdata.path.stem}.png"
plt.savefig(save_path, facecolor='black', dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="sum",
    table_name="filtered",
    cmap="plasma"
).pl.show(
    ax=ax, 
    title=f"Library Size per nucleus - {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/libsize_{sdata.path.stem}.png"
plt.savefig(save_path, facecolor='black', dpi=300, bbox_inches='tight')
plt.close()


fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="sum",
    table_name="filtered",
    cmap="plasma"
).pl.render_shapes(
    "intissue_poly", 
    outline=True, 
    outline_alpha=1, 
    fill_alpha=0,
    color="black",
).pl.show(
    ax=ax, 
    title=f"intissue exist - {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/try.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
#ax.set_facecolor('black')        # Set the plot area to black
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[5500, 11000],
    max_coordinate=[10000, 13000],
    target_coordinate_system=COORD_SYS,
).pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="labels_singler",
    table_name=TABLE_KEY,
).pl.show(
    ax=ax, 
    title=f"Annotation plot of {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/annotation_zoom_{sdata.path.stem}.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()

fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
#ax.set_facecolor('black')        # Set the plot area to black
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[10000, 4000],
    max_coordinate=[12000, 5500],
    target_coordinate_system=COORD_SYS,
).pl.render_shapes(
    SHAPES_KEY, 
    outline=True, 
    outline_alpha=0, 
    fill_alpha=1,
    color="labels_singler",
    table_name=TABLE_KEY,
).pl.show(
    ax=ax, 
    title=f"Annotation plot of {sdata.path.stem}", 
    coordinate_systems=COORD_SYS
)
save_path = f"{OUTPUT_DIR}/QC/annotation_zoom_{sdata.path.stem}.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close()
