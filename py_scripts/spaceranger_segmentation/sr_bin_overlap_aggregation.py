import spatialdata_plot
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import spatialdata as spd
from spatialdata_io import visium_hd

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_c26STAT3")

# i need the table of the bins to aggregate the nuclei, otherwise will not work.
# so i'll load everything with spd_io
wdata = visium_hd(path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs",
  dataset_id = "sr_bins",
  filtered_counts_file = False,
  bin_size = '002',
  bins_as_squares = True,
  annotate_table_by_labels = False,
  var_names_make_unique = True
)


s1_extent = sd.get_extent(sdata, coordinate_system='downscale_to_hires', exact=True, 
                  has_images=True, 
                  has_labels=True, 
                  has_points=True, 
                  has_shapes=True, 
                  elements=['blocco1_nuclei_intissue_boundaries']
)

def crop0(x,crs,bbox):
    return sd.bounding_box_query(
        x,
        min_coordinate=[bbox['x'][0], bbox['y'][0]],
        max_coordinate=[bbox['x'][1], bbox['y'][1]],
        axes=("x", "y"),
        target_coordinate_system=crs,
    )

wcrop = crop0(wdata, "sr_bins", s1_extent)

sdata['table_bins'] = wcrop['square_002um'].copy()

#' It's not working because I need to map bins table into the nuclei boundaries shapes... so funny!
#' I'm not sure if there is track of the couple cell_id -> location_id... 
#' For now I know only that some bin is shared between cell_ids



nuclei_key = "blocco1_c26STAT3_nuclei_boundaries"
"table_bins"
sdata["table_bins"] = sdata["table_bins"].reset_index()  # index becomes a column
# Filter bins intissue
bin_to_nuclei = sopa.spatial.sjoin(sdata,
  "table_bins",
  "blocco1_c26STAT3_nuclei_boundaries",
  how = "inner",
  predicate = "intersects",
  target_coordinate_system = "downscale_to_hires"
)

# Probably I need to fix the coordinate system of the new table... or not! or maybe I need to do it with geopandas

# there can be duplicated bins because inside 2 different tissue, I cannot exclude them all, i need to aggregate with sopa!
# bins_nuc = bins_nuc[bins_nuc['cell_id'].duplicated(keep=False) == False]
# bins_nuc = bins_nuc[['location_id', 'geometry', 'name']]

# Add in the sdata object both the bins and the intissue poly
sdata["bin_to_nuclei"] = bin_to_nuclei # ready to get inside the sdata 



# sopa metadata
sopa.utils.set_sopa_attrs(
  sdata,
  cell_segmentation_key = "blocco1_nuclei_hires_tissue_image",
  bins_table_key = "table_bins"
)

sopa.aggregate(sdata, key_added = 'nuclei_counts_nop', bins_key= "table_bins",
shapes_key = "blocco1_c26STAT3_nuclei_boundaries", expand_radius_ratio=0, min_transcripts=1,
min_intensity_ratio=0.2, no_overlap = True)


