import spatialdata as sd
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
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
from skimage.measure import regionprops_table
from py_scripts.segmentation.segm_functions import (assign_bins_to_nearest_nucleus, nuclei_filtering) 


# Read sdata with all the 'blocchi' (all in one or one sdata per blocco, not yet understood)

# segmentation and post processing

# output and clean
