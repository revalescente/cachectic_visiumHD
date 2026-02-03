import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import spatialdata_plot
import matplotlib.pyplot as plt
# from spatialdata.models import ShapesModel
# from spatialdata.transformations import Identity, get_transformation
# from skimage.measure import regionprops_table
# import sopa

sdata = sd.read_zarr("/mnt/europa/shared/sandri_visiumHD_data/bins/version_1.0.0/sdatas/blocco4_c26")

genes_faps = ["Pdgfra", "Ly6a", "Dcn"]
gene_trim = "Trim63"
gene_myh4 = "Myh4"

# to have invisible shapes where the values == 0 use this cmap
#user_def_cmap = "cividis"
user_def_cmap = "plasma"
new_cmap = set_zero_in_cmap_to_transparent(cmap=user_def_cmap)

shape_name = "intissue_008um"

# 1. Access the specific table from your sdata object
# Assuming you want to modify 'square_008um'
adata = sdata.tables['square_008um']

# 2. Replicate R: spe$faps <- colSums(counts(spe[genes_faps, ])) > 0
# We sum across the gene axis (axis=1) for the selected genes
faps_mask = adata[:, genes_faps].X.sum(axis=1).A1 > 0
adata.obs['faps'] = faps_mask

# 3. Replicate R: counts(spe)["Gene", ]
# Extract raw counts for specific genes
adata.obs['trim'] = adata[:, "Trim63"].X.toarray().flatten()
adata.obs['myh4'] = adata[:, "Myh4"].X.toarray().flatten()

# 4. Replicate R: df$faps_plot <- ifelse(df$faps, TRUE, NA)
# We use np.where to handle the conditional mapping
adata.obs['faps_plot'] = np.where(adata.obs['faps'], True, np.nan)

# 5. Replicate R: log1p transformations
# Note: log1p is natural log(1 + x)
adata.obs['trim_plot'] = np.log1p(adata.obs['trim'])
adata.obs['myh4_plot'] = np.log1p(adata.obs['myh4'])

# 6. Update the table back in the sdata object (if not modified in-place)
sdata.tables['square_008um'] = adata

plt.figure(figsize=(20, 20))
ax = plt.gca()
# 1. Specify the column name in 'color'
# 2. Ensure 'groups' or 'palette' is used if it's categorical
sdata.pl.render_shapes(
    shape_name, 
    color="trim_plot",  # Add your column name here,
    table_name = "square_008um",
    outline=False, 
    #outline_alpha=1, 
    #outline_width=0, 
    fill_alpha=1,               # Increased slightly to see the 'new column' data
    cmap=new_cmap
).pl.show(
    ax=ax, 
    coordinate_systems="blocco4", 
    save="output_python/grafico_emma.png"
)
 
fig = plt.figure(figsize=(20, 20))
ax = plt.gca()
ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    shape_name, 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="trim_plot",
    table_name="square_008um",
    cmap=new_cmap
).pl.show(
    ax=ax, 
    title=f"Aspect Ratio {sdata.path.stem}", 
    coordinate_systems="blocco4"
)
save_path = "/mnt/europa/valerio/figures/grafico_emma2.png"
plt.savefig(save_path, facecolor='black', dpi=300, bbox_inches='tight')
plt.close()
