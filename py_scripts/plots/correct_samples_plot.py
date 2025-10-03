import spatialdata as sd
import spatialdata_plot
from spatialdata import SpatialData
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import sopa
import json
from shapely.affinity import scale
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import py_scripts.pp_sdata.pp_functions as pp


# dictionary to manage the various samples
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    blocco_sample_bbox_dict = json.load(f)
    
for blocco, samples in blocco_sample_bbox_dict.items():
  for sample,_ in samples.items():
    sdata = read_zarr_standardized(f"/mnt/europa/valerio/data/zarr_store/blocchi/{blocco}_{sample}.zarr")
    plt.figure(figsize=(20, 20))
    ax = plt.gca()
    sdata.pl.render_images(f"{blocco}_full_image", scale = "scale2").pl.render_shapes(f"{blocco}_intissue", outline=True, outline_alpha=1, outline_width=3, fill_alpha=0
    ).pl.show(ax = ax, coordinate_systems = blocco, save = f"output_python/testing/{blocco}_{sample}.png")
 

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco4_c26.zarr")
sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco9_c26SMAD23.zarr")

plt.figure(figsize=(50, 50))
ax = plt.gca()
# sdata.query.bounding_box(
#     axes=["x", "y"],
#     min_coordinate=[5500, 14500],
#     max_coordinate=[6000, 15000],
#     target_coordinate_system="blocco9_c26SMAD23",
# )
sdata.pl.render_images("blocco9_c26SMAD23_full_image", scale = "scale2"
).pl.render_shapes("blocco9_c26SMAD23_filtered_nuclei", outline=True, outline_alpha=1, outline_width=3, fill_alpha=0
).pl.show(ax = ax, coordinate_systems = "blocco9_c26SMAD23", save = "output_python/presentation/b9_c26SMAD23_result.png")






















