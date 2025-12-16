import matplotlib.pyplot as plt
import numpy as np
from skimage import io

# 1. Load the image
path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/blocco2/outs/spatial/tissue_hires_image.png"
img = io.imread(path)

# 2. Downsample for faster plotting (Optional but HIGHLY recommended for 6k images)
# Taking every 10th pixel is enough to see the orientation.
# Remove [::10, ::10] if you really need full resolution.
preview = img[::10, ::10] 

# 3. Define the transformations
# Note: np.rot90 rotates Counter-Clockwise
transformations = [
    ("Original", preview),
    ("Rotated 90 CCW", np.rot90(preview, k=1)),
    ("Rotated 180", np.rot90(preview, k=2)),
    ("Rotated 270 CCW (-90)", np.rot90(preview, k=3)),
    ("Flip Left-Right (Mirror)", np.fliplr(preview)),
    ("Flip Up-Down", np.flipud(preview)),
    ("Flip LR + Rot 90", np.rot90(np.fliplr(preview), k=1)), # Common in microscopy
    ("Transpose (Swap X/Y)", np.transpose(preview, (1, 0, 2)))
]

# 4. Create a grid plot (2 rows, 4 columns)
fig, axes = plt.subplots(2, 4, figsize=(24, 12))
axes = axes.flatten()

for ax, (name, transformed_img) in zip(axes, transformations):
    ax.imshow(transformed_img)
    ax.set_title(name, fontsize=14)
    ax.axis('off') # Hide axes for cleanliness

# 5. Save and close (Do not show)
output_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/roba_da_buttare/orientation_check.png"
plt.tight_layout()
plt.savefig(output_path, dpi=150)
plt.close(fig) # Frees up memory immediately

print(f"Saved comparison grid to: {output_path}")

# ------------------------------------------------------------------------------

# so let's do it with the spatialdata object, insert the geojson with an affine transformation: 180 degrees rotation

import spatialdata as sd
import spatialdata_plot
#from spatialdata import SpatialData
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import (Identity, set_transformation, Scale, Affine)
import matplotlib.pyplot as plt
#import squidpy as sq
import numpy as np
#import scanpy as sc
import geopandas as gpd
import pandas as pd
#import os
import sopa
#from pathlib import Path
#from sopa.io.standardize import sanity_check
#from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
#from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf
from py_scripts.utils.utils_fun import read_from_json

blocco = 'blocco2'
samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]
# data
sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[2]}_storto.zarr")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[1]}_storto.zarr")
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[0]}_storto.zarr")

geojson_path = f"/mnt/europa/valerio/data/json/geojson_dir/GFP_inflamed/{blocco}_GFP_rotated.geojson"
gfp_poly = gpd.read_file(geojson_path)
gfp_poly = gfp_poly.set_crs(None, allow_override=True)

# devo inserire la gfp_poly rotata di 180 gradi, quindi devo definire una trasformazione affine di quel tipo
# 
# theta = math.pi / 6
# rotation = Affine(
#     [
#         [math.cos(theta), -math.sin(theta), 0],
#         [math.sin(theta), math.cos(theta), 0],
#         [0, 0, 1],
#     ],
#     input_axes=("x", "y"),
#     output_axes=("x", "y"),
# )
# 
# set_transformation(sdata.images["raw_image"], rotation, to_coordinate_system="global")

#gfp_parse = ShapesModel.parse(gfp_poly, transformations={blocco: Identity()})


# sdata2.shapes["GFP_poly"] = gfp_parse
# sdata3.shapes["GFP_poly"] = gfp_parse

angle = np.pi  # 30 degrees in radians
cos_angle = np.cos(angle)
sin_angle = np.sin(angle)
affine = Affine(
    [[cos_angle, -sin_angle, 0], [sin_angle, cos_angle, 0], [0, 0, 1]], input_axes=("x", "y"), output_axes=("x", "y")
)

# Step 3: Put these squares into a GeoPandas GeoDataFrame
gfp_parse = ShapesModel.parse(gfp_poly, transformations={"blocco2": affine})
gfp_parse 

gfp_sdata1 = gfp_parse[gfp_parse['name'].str.contains("infiammazione_gr1", case=False, na=False)]

sdata1.shapes["GFP_poly"] = gfp_parse
sdata1.shapes["GFP_poly"] = gfp_sdata1

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_shapes('GFP_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_gfp_affine_trans.png')


# retry

gfp_parse = ShapesModel.parse(gfp_poly, transformations={"global": Identity()})
sdata1.shapes["GFP_poly"] = gfp_parse
sdata2.shapes["GFP_poly"] = gfp_parse
sdata3.shapes["GFP_poly"] = gfp_parse

theta = np.pi
rotation = Affine(
    [
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ],
    input_axes=("x", "y"),
    output_axes=("x", "y"),
)

set_transformation(sdata1.shapes["GFP_poly"], rotation, to_coordinate_system="blocco2")
set_transformation(sdata2.shapes["GFP_poly"], rotation, to_coordinate_system="blocco2")
set_transformation(sdata3.shapes["GFP_poly"], rotation, to_coordinate_system="blocco2")

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_shapes('GFP_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_gfp_affine_trans.png')

# non va neanche così... testiamo con intissue_poly la rotazione

set_transformation(sdata1.shapes["intissue_poly"], rotation, to_coordinate_system="blocco2")

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('intissue_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=1
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = f'output_python/roba_da_buttare/prova_rotazione.png'
)

#' no la rotazione lo porta completamente fuori, devo farla rispetto al centro dell'immagine però.... 
#' ma io ho l'immagine del blocco intero dio cane... 
#' 

# ok allora farò cosi, inserisco l'immagine gfp_poly in un sdata tutto nuovo insieme all'intissue_poly originale, 
# così faccio combaciare loro, o almeno ci provo e sono più leggeri da ritagliare e manovrare

#' ANZI, prima provo a girare l'immagine del blocco2 e allinearla così alle shapes, 
#' forse è più facile girare l'immagine che i poligoni... dio cane perché non c'ho pensato prima porco dio

sdata1['full_image']

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('intissue_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=1
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = f'output_python/roba_da_buttare/prova_rotazione.png'
)


# lastly, rotate the gpd from the center... and hope
gfp_co = gfp_poly.copy()
gfp_co['geometry'] = gfp_co.rotate(180, origin='center')

gfp_co_parse = ShapesModel.parse(gfp_co, transformations={"blocco2": Identity()})
gfp_co_parse 

sdata1.shapes["GFP_poly"] = gfp_co_parse
sdata2.shapes["GFP_poly"] = gfp_co_parse
sdata3.shapes["GFP_poly"] = gfp_co_parse
