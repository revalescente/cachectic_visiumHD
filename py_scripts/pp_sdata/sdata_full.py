import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
import matplotlib.pyplot as plt
import squidpy as sq
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
import sopa

# reading data with visium reader
# spe_blocco1 = visium_hd(path = '/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/blocco1/outs', 
#                 dataset_id='blocco1', 
#                 filtered_counts_file=False, 
#                 bin_size='002', 
#                 bins_as_squares=True, 
#                 annotate_table_by_labels=False, 
#                 fullres_image_file='/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco1_20x.tif', 
#                 load_all_images=False, 
#                 var_names_make_unique=True)

block_numbers = [1, 2, 3, 4, 5, 6, 7, 9]
spe_blocks = {}

for i in block_numbers:
    block_name = f"spe_blocco{i}"
    spe_blocks[block_name] = visium_hd(
        path=f"/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/blocco{i}/outs",
        dataset_id=f"blocco{i}",
        filtered_counts_file=False,
        bin_size='002',
        bins_as_squares=True,
        annotate_table_by_labels=False,
        fullres_image_file=f"/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco{i}_20x.tif",
        load_all_images=False,
        var_names_make_unique=True
    )

# Assuming spe_blocks is your dictionary of datasets
for block_name, spe in spe_blocks.items():
    spe.write(f"/mnt/europa/valerio/data/{block_name}.zarr")

# preprocess function

def sdata_pp(scale_factor=3.78633328):
  
  # read sdata
  spe = sd.read_zarr('/mnt/europa/valerio/data/spe_blocco1.zarr')
  
  # read the intissue polygons
  intissue_poly = gpd.read_file('~/geojson_dir/tissue_hires_image_blocco1.geojson')
  intissue_scaled = intissue_poly.set_crs(None, allow_override=True).copy()

  # Apply scaling to all geometries in intissue_poly
  xfact = scale_factor
  yfact = scale_factor
  intissue_scaled['geometry'] = intissue_poly['geometry'].apply(
      lambda geom: scale(geom, xfact=xfact, yfact=yfact, origin=(0, 0))
  )
  
  # extract bins shapes keeping the index 'location_id' for the filtering
  bins = spe['blocco1_square_002um'].reset_index()  # location_id becomes a column
  # filter bins intissue
  intersection = gpd.overlay(bins, intissue_scaled, how='intersection')
  # add in the spe object
  intersection_parse = ShapesModel.parse(intersection, transformations={'blocco1': Identity()})
  spe.shapes['intissue_002um'] = intersection_parse


# read and print the Zarr data
spe = sd.read_zarr('spe_blocco1_mod001.zarr')
print(spe)



















