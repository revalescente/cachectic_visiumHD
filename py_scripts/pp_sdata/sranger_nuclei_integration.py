import re
import json
from shapely.affinity import scale
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import geopandas.testing 
import sopa
from sopa.io.standardize import sanity_check
import spatialdata as sd
import spatialdata_plot
from spatialdata.transformations import Identity
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from py_scripts.utils.utils_fun import read_from_json
import scanpy as sc
import os

# file needed for the single sample prep
sdata = sd.read_zarr('/mnt/europa/valerio/data/zarr_store/samples/blocco1_sham.zarr')
geojson_path = '/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/nucleus_segmentations.geojson'
parquet_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/barcode_mappings.parquet"
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')
nuclei_matrix_path = '/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/raw_feature_cell_matrix'

# nuclei boundaries integrated
sranger_nuclei = geojson_preparation(geojson_path = geojson_path, samples_dict = samples_dict, sample_key = 'c26foxO', blocco_key ='blocco1')
sdata['sranger_nuclei'] = sranger_nuclei

# adata intergation

sdata['sranger_table'] = table_integration(sranger_nuclei_gpd = sranger_nuclei, nuclei_matrix_path = nuclei_matrix_path, region_key = 'sranger_nuclei')



