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
import shutil
from shapely.affinity import scale
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent

# leggiamo come lista
samples = ["blocco1_c26foxO", "blocco1_sham", "blocco3_c26murf1"]
spe_blocks_list = []

for sample in samples:
    spe_blocks_list.append(sd.read_zarr(f'/mnt/europa/valerio/data/zarr_store/blocchi/{sample}.zarr'))

# leggiamo come dizionario 
block_numbers = [4, 7, 9]
spe_blocks = {}
for i in block_numbers:
  block_name = f"blocco{i}"
  spe_blocks[block_name] = sd.read_zarr(f'/mnt/europa/valerio/data/spe_blocco{i}.zarr')

# concateniamo i due e vediamo le differenze
spe_list = sd.concatenate(spe_blocks_list, concatenate_tables = True)
spe_dict = sd.concatenate(spe_blocks, concatenate_tables = True)
