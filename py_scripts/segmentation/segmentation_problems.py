from scipy.ndimage import gaussian_filter
import spatialdata as sd
#import spatialdata_plot
#from spatialdata import SpatialData
#from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
#from spatialdata.transformations import (Identity, set_transformation)
import numpy as np
#import scanpy as sc
#import matplotlib.pyplot as plt
import geopandas as gpd
import pandas as pd
import os
#import gc
import sopa
import json
#import time
#from pathlib import Path
from sopa.io.standardize import sanity_check, read_zarr_standardized
#from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
#from skimage.measure import regionprops_table
import py_scripts.segmentation.segm_functions as sf

sopa.settings.auto_save_on_disk = False

# Probabilmente questo campione non funziona perché c'è almeno una patches con una dimensione
# pari a zero (altezza o larghezza). Gustave sta investigando.

sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco6_sham.zarr")

spatial_image = sopa.get_spatial_image(sdata)
image_channels = spatial_image.coords["c"].values
image_shapes_vector = []

for patch in sdata["image_patches"].geometry.values:
  bounds = [int(x) for x in patch.bounds]
  image = spatial_image.sel(
        c=image_channels,
        x=slice(bounds[0], bounds[2]),
        y=slice(bounds[1], bounds[3]),
  ).values
  image = np.stack([gaussian_filter(c, sigma=1) for c in image])
  image_shapes_vector.append(image.shape)
  
df = pd.DataFrame(image_shapes_vector, columns=['channels', 'height', 'width'])

df.describe()
np.sum(df['height'] < 80 )
np.sum(df['height'] < 40 )
np.sum(df['height'] < 20 )
np.sum(df['width'] < 80 )
np.sum(df['width'] < 40 )
np.sum(df['width'] < 20 )
np.sum(df['width'] == 00 )


# -----------------------------------------------------------------------------

# NO overlap - expand ratio != 0 non funziona - estraiamo i parquet
# per investigare la forma dei nuclei. 

sdata = read_zarr_standardized("/mnt/europa/valerio/data/zarr_store/blocchi/blocco2_c26.zarr")

sopa.aggregate(sdata, key_added = 'nuclei_counts_exp1', bins_key= "filtered",
  shapes_key = "blocco2_nuclei_boundaries", 
  expand_radius_ratio=1, 
  min_transcripts=1,
  min_intensity_ratio=0.1, 
  no_overlap = True
)







sopa.settings.auto_save_on_disk = True
