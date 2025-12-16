import spatialdata as sd
import spatialdata_plot
#from spatialdata import SpatialData
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import (Identity, set_transformation, Scale)
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
sdata = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[2]}_storto.zarr")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[1]}_storto.zarr")
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/old_samples/{sample_keys[0]}_storto.zarr")
sdata_prova = sdata.copy()
del sdata_prova['nuclei_boundaries']
sdata_prova = sf.segmentation_step(sdata_prova)

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('full_image', scale = "scale3"
).pl.render_shapes('nuclei_boundaries'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_segments.png'
)

del sdata_prova['nuclei_counts_nop']
del sdata_prova['filtered_nuclei']
sdata_prova = sf.postprocess_step(sdata_prova, expand_radius_ratio = 0, min_transcripts = 4, min_intensity_ratio = 0.15, no_overlap = True)

# provo a riaggiungere l'immagine e vediamo se cosi sta dritta?? bo

# no perché torna tutto l'unico problema è aggiungere la GFP e le aree da rimuovere... dio cane
geojson_path = f"/mnt/europa/valerio/data/json/geojson_dir/GFP_inflamed/{blocco}_GFP_and_inflamed.geojson"
gfp_poly = gpd.read_file(geojson_path)
gfp_poly = gfp_poly.set_crs(None, allow_override=True)
gfp_parse = ShapesModel.parse(gfp_poly, transformations={blocco: Identity()})

sdata.shapes["GFP_poly"] = gfp_parse

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('GFP_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems=blocco, 
    save = f'output_python/roba_da_buttare/{blocco}_gfp_poly.png'
)
gr1_poly = gfp_poly[gfp_poly['name'].str.contains("gr1", case=False, na=False)]
# --------------------------
# roba di emma: vediamo se la separazione dei campioni e' avvenuta correttamente ecccc, non centra con cio che c'e sopra.

blocco_key = 'blocco1'
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]
# data
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/binned_samples/{sample_keys[2]}")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/binned_samples/{sample_keys[1]}")
sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/binned_samples/{sample_keys[0]}")

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_shapes('nuclei_boundaries'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_segments.png'
)

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata2.pl.render_shapes('nuclei_boundaries'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_segments.png'
)

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata3.pl.render_shapes('nuclei_boundaries'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_segments.png'
)


# blocco 2 capiamo perché l'immagine è capovolta e perché è così capovolta

from skimage import io

img = io.imread("/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/blocco2/outs/spatial/tissue_hires_image.png")

# Transformation for shapes (from pixel to downscale_to_hires)
image_transform = {
    "blocco2": Scale(np.array([3.56899999, 3.56899999]), axes=("x", "y"))
}
img_parsed = Image2DModel.parse(img, dims = ["y", "x", "c"], transformations=image_transform)

sdata["hires_due"] = img_parsed

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('hires_due'
).pl.render_shapes('intissue_008um'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_prova_img_dritta.png'
)

# è storta rispetto alle geometrie

# try full res and see
fullimg = io.imread("/mnt/europa/valerio/HE_images/color_corrected/pp_blocco2_20x.tif")
fullimg_parsed = Image2DModel.parse(fullimg, dims = ["y", "x", "c"], transformations={"blocco2" : Identity()}, scale_factors = [2, 2, 2])
sdata["fullres"] = fullimg_parsed

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('fullres', scale = "scale3"
).pl.render_shapes('intissue_008um'
).pl.show(ax = ax, coordinate_systems="blocco2", 
    save = 'output_python/roba_da_buttare/blocco2_prova_fullimg_dritta.png'
)


# -----------------------------

plt.figure(figsize=(50, 50))
ax = plt.gca()
temp_sd.pl.render_images('fluo_image', cmap = 'grey', scale = 'scale3'
).pl.show(ax = ax, coordinate_systems='blocco2', save = f'output_python/prova_gfp_image.png')


plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images('full_image', scale = 'scale3'
).pl.render_images('fluo_image', cmap = 'grey', scale = 'scale3'
).pl.show(ax = ax, coordinate_systems='blocco2', save = f'output_python/prova_gfp_image.png')



