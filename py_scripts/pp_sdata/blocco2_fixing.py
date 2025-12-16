import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
from spatialdata import SpatialData
import matplotlib.pyplot as plt
import numpy as np
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

# per il blocco2 non devo scalare gli intissue

blocco2 = visium_hd(
        path="/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/blocco2/outs",
        dataset_id="blocco2",
        filtered_counts_file=False,
        bin_size='002',
        bins_as_squares=True,
        annotate_table_by_labels=False,
        fullres_image_file="/mnt/europa/valerio/HE_images/color_corrected/pp_blocco2_20x.tif",
        load_all_images=False,
        var_names_make_unique=True
    )


# plt.figure(figsize=(20, 20))
# ax = plt.gca()
# blocco2.pl.render_images('blocco2_full_image', scale = 'scale3').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_control.png')

geojson_dir = "/mnt/europa/valerio/data/json/geojson_dir/"
blocco_key = 'blocco2'
hires_key = 'blocco2_hires_image'
geojson_path = f"{geojson_dir}tissue_fullres_image_{blocco_key}.geojson"

# Read intissue polygons
intissue_poly = gpd.read_file(geojson_path)
intissue_poly = intissue_poly.set_crs(None, allow_override=True)

# Add intissue poly in the sdata
intissue_rename = "intissue_poly"
intissue_parse = ShapesModel.parse(intissue_poly, transformations = {blocco_key: Identity()})
blocco2.shapes[intissue_rename] = intissue_parse
    
# Extract bins shapes keeping the index 'location_id' for the filtering
bins_key = [key for key in blocco2.shapes if re.search(r'_square_002um$', key)][0]
blocco2[bins_key] = blocco2[bins_key].reset_index()  # location_id becomes a column
# Filter bins intissue
bins_intissue = sopa.spatial.sjoin(blocco2,
      bins_key,
      intissue_rename,
      how = "inner",
      predicate = "intersects",
      target_coordinate_system = blocco_key
)
    # there can be duplicated bins because inside 2 different tissue, I exclude them
bins_intissue = bins_intissue[bins_intissue['location_id'].duplicated(keep=False) == False]
bins_intissue = bins_intissue[['location_id', 'geometry', 'name']]
# Add in the sdata object both the bins and the intissue poly
bins_shape_rename = "intissue_002um"
blocco2.shapes[bins_shape_rename] = bins_intissue # ready to get inside the sdata 


# Annotate the table with the filtered bins
blocco2["square_002um"].obs["region"] = pd.Categorical([bins_shape_rename] * len(blocco2["square_002um"]))
blocco2["square_002um"].uns["spatialdata_attrs"] = {
        "region": bins_shape_rename,  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
}
blocco2.set_table_annotates_spatialelement("square_002um", region=bins_shape_rename)    
    
# Filter the table
blocco2['filtered'] = sd.match_table_to_element(blocco2, 
        element_name=bins_shape_rename, 
        table_name="square_002um"
)

# Map the 'exp_condition' in the 'filtered' table
location_to_name = blocco2[bins_shape_rename].set_index('location_id')['name']
# Map the names to adata.obs using location_id
blocco2['filtered'].obs['sample_id'] = blocco2['filtered'].obs['location_id'].map(location_to_name)

# create a new sdata with only the interesting elements.
final = blocco2.subset([bins_shape_rename, intissue_rename, 'filtered'], filter_tables = False)
image_rename = "full_image"
final[image_rename] = blocco2[f'{blocco_key}_full_image'].copy()
# 
# plt.figure(figsize=(20, 20))
# ax = plt.gca()
# final.pl.render_images('full_image', scale = 'scale3').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_control.png')

# b2_old = sd.read_zarr('/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco2.zarr')
# prima di salvare controlliamo che sia tutto ok e che i due oggetti conservino le stesse informazioni. che non sia mai.
# ok sono uguali. l'immagine in final e' dritta quindi salvo e procedo a dividere i campioni.
final.write('/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco2.zarr')

# ------

# divisione campioni
samples_dict = pp.read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')
subset_dict = {name: samples_dict[name] for name in ['blocco2'] if name in samples_dict}

# blocco2 = sd.read_zarr('/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco2.zarr')

plt.figure(figsize=(20, 20))
ax = plt.gca()
final.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=[0,13000],
                max_coordinate=[21414, 20069],
                target_coordinate_system='blocco2'
            ).pl.render_images('full_image', scale = 'scale3').pl.render_shapes('intissue_poly').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_control_c26.png')


input_dir  = "/mnt/europa/valerio/data/zarr_store/filtered/"
output_dir = "/mnt/europa/valerio/data/zarr_store/"

    # "blocco2": {
    #     "c26": {
    #         "min_coordinate": [0,13000],
    #         "max_coordinate": [21414, 20069],
    #         "sample_key": "blocco2_c26"
    #     },
    #     "c26murf1": {
    #         "min_coordinate": [0, 8750],
    #         "max_coordinate": [21414, 13750],
    #         "sample_key": "blocco2_c26murf1"
    #     },
    #     "c26SMAD23": {
    #         "min_coordinate": [0, 0],
    #         "max_coordinate": [21414, 8750],
    #         "sample_key": "blocco2_c26SMAD23"
    #     }

sdata = sd.read_zarr(f"{input_dir}filtered_blocco2.zarr")

# Ho fatto un errore sui poligony intissue, ho chiamato un tessuto nel modo sbagliato. non e' foxO
# sdata['intissue_poly'].loc[sdata['intissue_poly']['name'] == 'c26foxO', 'name'] = 'c26murf1'
# print(sdata['intissue_poly']['name'])

sample = 'c26'
min_coordinate = [0,13000]
max_coordinate = [21414, 20069]
# Bounding box query
sdata_bbox = sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coordinate,
                max_coordinate=max_coordinate,
                target_coordinate_system="blocco2"
)

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata_bbox.pl.render_images('full_image', scale = 'scale3').pl.render_shapes('intissue_poly').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_control_murf1_v2.png')


# --- INTEGRATION: filter table bins and intissue_poly ---
# let's use get_value and filter the table only, then match both the elements
filty = sd.get_values(
              value_key='sample_id', 
              element=sdata_bbox['filtered']
)['sample_id']
sdata_bbox["filtered"] = sdata_bbox["filtered"][filty == sample]
            
# match elements
sdata_bbox['intissue_002um'] = sd.match_element_to_table(
              sdata_bbox, 
              element_name='intissue_002um', 
              table_name='filtered'
)[0]['intissue_002um']

# filtering 'intissue' by exp_condition (column "name") 
sdata_bbox['intissue_poly'] = sdata_bbox['intissue_poly'][sdata_bbox['intissue_poly']['name'] == sample]

# setup sopa metadata 
sopa.utils.set_sopa_attrs(
                sdata_bbox,
                cell_segmentation_key = "full_image",
                tissue_segmentation_key = "full_image",
                bins_table_key = "filtered"
)
            

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata_bbox.pl.render_images('full_image', scale = 'scale3').pl.render_shapes('intissue_poly').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_control_c26_v3.png')


# Sanity check
try:
  sanity = sanity_check(sdata_bbox)
  if sanity is None:
    out_path = f"{output_dir}blocco2_{sample}.zarr"
    sdata_bbox.write(out_path)
    print(f"Wrote {out_path}")
  else:
     print(f"Sanity check returned a value for {blocco}_{sample}; skipping write.")
except AssertionError as e:
  print(f"Sanity check failed for {blocco}_{sample}: {e}")


# ------------------------------------------------------------------------------

# inserire i nuclei segmentati e la tabella? o devo rifare tutto?cazzzzzoooooooo
input_dir  = "/mnt/europa/valerio/data/zarr_store/samples/"
output_dir = "/mnt/europa/valerio/data/zarr_store/"

sdata1 = sd.read_zarr(f"{input_dir}blocco2_c26.zarr") # sbagliato (full image storta non si sa come)
sdata2 = sd.read_zarr(f"{output_dir}blocco2_c26.zarr") # nuovo 

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_images('full_image', scale = 'scale3').pl.render_shapes('filtered_nuclei').pl.render_shapes('intissue_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/roba_da_buttare/b2_c26_control_vecchio.png')

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata2.pl.render_images('full_image', scale = 'scale3').pl.render_shapes('intissue_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems="blocco2_c26", save = 'output_python/roba_da_buttare/b2_c26_control_nuovo.png')

