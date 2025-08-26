from spatialdata_io import visium_hd
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.models._utils import get_channel_names, set_channel_names
from spatialdata import SpatialData
import spatialdata as sd
from spatialdata.transformations import (
    Identity,
    get_transformation,
    set_transformation,
)
from shapely.affinity import scale
from shapely import union_all
from anndata import AnnData
from skimage import imread
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
import math

# to add elements in the zarr store:
# sdata.write_element(["element1", "element2"])
# delete element from memory (from sdata open in python)
# del sdata['element1']
# delet element from disk
# sdata.delete_element_from_disk("element_name")

# reading data with visium reader
spe = visium_hd(path = '/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/blocco1/outs', 
                dataset_id='blocco1', 
                filtered_counts_file=False, 
                bin_size='002', 
                bins_as_squares=True, 
                annotate_table_by_labels=False, 
                fullres_image_file='/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE/Project_Blocco1_20x.tif', 
                load_all_images=False, 
                var_names_make_unique=True)

# write the data to disk
spe.write('spe_blocco1.zarr')

# read and print the Zarr data
spe = sd.read_zarr('spe_blocco1_mod001.zarr')
print(spe)

# read the intissue polygons
intissue_poly = gpd.read_file('~/data/geojson_dir/tissue_hires_image_blocco1.geojson')
intissue_poly = intissue_poly.set_crs(None, allow_override=True)

# Apply scaling to all geometries in intissue_poly
xfact =  3.78633328  
yfact =  3.78633328

intissue_scaled = intissue_poly.copy()
intissue_scaled['geometry'] = intissue_scaled['geometry'].apply(
    lambda geom: scale(geom, xfact=xfact, yfact=yfact, origin=(0, 0))
)

# extract bins shapes keeping the index 'location_id' for the filtering
bins = spe['blocco1_square_002um'].reset_index()  # location_id becomes a column
# filter bins intissue
intersection = gpd.overlay(bins, intissue_scaled, how='intersection')
# add in the spe object
intersection_parse = ShapesModel.parse(intersection, transformations={'blocco1': Identity()})
spe.shapes['intissue_002um'] = intersection_parse

# add intissue_002um to the coord system 'blocco1'
transformed = get_transformation(spe.images['blocco1_full_image'], to_coordinate_system='blocco1')
set_transformation(spe.shapes['intissue_002um'], transformed, to_coordinate_system='blocco1')

# filter genes counts table by row in tissue
counts = spe.tables['square_002um']

print(len(set(counts.obs.index) & set(intersection_parse['location_id'].astype(str).unique())))

location_ids = intersection_parse['location_id'].astype(str).unique()
filtered_counts = counts[counts.obs['location_id'].astype(str).isin(location_ids)].copy()
filtered_counts.obs['exp_cond'] = intersection_parse['name'].values

# connect the AnnData with the Shapes object
counts_parse = TableModel.parse(filtered_counts)

# setting up metadata

# spe['square_002um'].uns :
# {'spatialdata_attrs': {'instance_key': 'location_id', 'region': 'blocco1_square_002um', 'region_key': 'region'}}

counts_parse.uns["spatialdata_attrs"] = {
    "region": "intissue_002um",  # name of the Shapes element we will use later 
    "region_key": "region",  # column in adata.obs that will link a given obs to the elements it annotates
    "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
}
# all the rows of adata annotate the same element, called "intissue_002um" (as we declared above)
counts_parse.obs["region"] = pd.Categorical(["intissue_002um"] * len(counts_parse))
counts_parse.obs[["region", "location_id"]]

spe.tables['filtered'] = counts_parse

# subsetting spe to obtain the filtered sdata with fewer objects
filtered = spe.subset(['blocco1_full_image', 'blocco1_hires_image', 'blocco1_intissue', 'intissue_002um', 'filtered'])
filtered

# adjusting coord systems
transf = get_transformation(spe.images['blocco1_hires_image'], to_coordinate_system="blocco1")
# Set a transformation for an element to a coordinate system
set_transformation(filtered.shapes['blocco1_intissue'], transf, to_coordinate_system="blocco1")

# save the new object
filtered.write("b1_filtered.zarr")


# Read filtered and cleaned data
filtered = sd.read_zarr(os.path.expanduser("~/data/b1_cleaned.zarr")) # non mi chiedere perché 

# add convolved image
conv_img = imread(os.path.expanduser("~/HE_images/preprocessed/convolved_blocco1_20x.tif"))
conv_img = np.expand_dims(image, axis=0)
conv_img_parsed = Image2DModel.parse(conv_img, dims=("c", "y", "x"),  transformations={'blocco1': Identity()}, c_coords = "conv")
filtered['blocco1_conv_image'] = conv_img_parsed
filtered.write_element("blocco1_conv_image")
