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
import shutil
from shapely.affinity import scale
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from collections import Counter

# dizionario campioni
with open('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json', 'r') as f:
    blocco_sample_bbox_dict = json.load(f)

# extract the sample_key in a list
samples = [
    details['sample_key'] 
    for blocco in blocco_sample_bbox_dict.values() 
    for details in blocco.values()
]

# leggiamo come lista
# samples = ["blocco1_c26foxO", "blocco1_sham", "blocco3_c26murf1"]
spe_blocks_list = []

for sample in samples:
    spe_blocks_list.append(sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/samples/{sample}.zarr'))

# leggiamo come dizionario 
spe_blocks = {}
for sample in samples:
  spe_blocks[sample] = sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/blocchi/{sample}.zarr')

# i raster nuclei devo eliminarli
for el,sample in zip(raster_nuclei_vector,spe_blocks_list):
  del sample[el]

# concateniamo i due e vediamo le differenze
spe_list = sd.concatenate(spe_blocks_list, concatenate_tables = True)
spe_dict = sd.concatenate(spe_blocks, concatenate_tables = True)

spe_list.write('/mnt/europa/valerio/data/zarr_store/concat_all_samples_almost.zarr')
















