import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
from spatialdata import SpatialData
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import sopa
import re
import json
from shapely.affinity import scale
from sopa.io.standardize import sanity_check, read_zarr_standardized
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
 


























