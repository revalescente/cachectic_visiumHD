# Let's try to extract features from the nuclei segmentation geometries (or labels)
# Since I don't know how to make it work with spatialdata

import spatialdata_plot
import matplotlib.pyplot as plt
import squidpy as sq
import os
import sopa
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized

b1_stat3 = read_zarr_standardized(os.path.expanduser("~/data/b1_stat3.zarr"))

sopa.utils.set_sopa_attrs(
    b1_stat3,
    cell_segmentation_key="blocco1_full_image",
    tissue_segmentation_key="blocco1_full_image",
    bins_table_key="cleaned",
    boundaries_key = "blocco1_intissue"
)


