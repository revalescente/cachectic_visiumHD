import re
import os
import json
from shapely.affinity import scale
import pandas as pd
import geopandas as gpd
import sopa
from sopa.io.standardize import sanity_check, read_zarr_standardized
import spatialdata as sd
from spatialdata_io import visium_hd
from spatialdata.transformations import Identity
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from py_scripts.utils.utils_fun import read_from_json

def preprocess_step(
    all=True, 
    block_numbers=None, 
    input_path="/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/",
    output_path="/mnt/europa/valerio/data/zarr_store/",
    images_path="/mnt/europa/valerio/HE_images/color_corrected/"
):
    """
    Function to prepare data for segmentation: filtering useless things.

    If you don't specify which blocks to process, the function will process all of them.
    If all=True, runs both reading and filtering steps.
    If all=False, runs only the filtering step.
    """
    # Which blocks to preprocess
    if block_numbers is None:
        block_numbers = [1,2,3,4,5,6,7,9]

    for i in block_numbers:
        blocco_key = f"blocco{i}"

        # 1. Reading data from source (only if all=True)
        if all:
            try:
                sdata = visium_hd(
                    path=f"{input_path}blocco{i}/outs",
                    dataset_id=blocco_key,
                    filtered_counts_file=False,
                    bin_size='002',
                    bins_as_squares=True,
                    annotate_table_by_labels=False,
                    fullres_image_file=f"{images_path}pp_blocco{i}_20x.tif",
                    load_all_images=False,
                    var_names_make_unique=True
                )
                sdata.write(f"{output_path}general/{blocco_key}.zarr")
                del sdata
            except Exception as e:
                print(f"Error reading {blocco_key}: {e}")

        # 2. Filtering of blocks (always runs)
        try:
            result = sdata_pp(f"{output_path}general/{blocco_key}.zarr")
            print(f"Processed {blocco_key}")
            # Resaving as filtered data
            result.write(f'{output_path}filtered/filtered_{blocco_key}.zarr')
            print(f"Zarr Store created for {blocco_key}")
        except Exception as e:
            print(f"Error processing {blocco_key}: {e}")
    print("Job done")

# ------------------------------------------------------------------------------

def sdata_pp(path = None,
             geojson_dir = "/mnt/europa/valerio/data/json/geojson_dir/"
):
    # Infer blocco name from path
    # e.g. '/mnt/europa/valerio/data/zarr_store/general/blocco1.zarr' -> 'blocco1'
    blocco_match = re.search(r'(blocco\d+)', path)
    if not blocco_match:
        raise ValueError("Cannot infer blocco name from path.")
    blocco_key = blocco_match.group(1)

    # Read spatialdata object
    sdata = sd.read_zarr(path)

    # Find hires image key using regex
    hires_keys = [key for key in sdata.images if re.search(r'_hires_image$', key)]
    if not hires_keys:
        raise ValueError("No hires image found in sdata.images")
    hires_key = hires_keys[0]

    # Access the y scale value (index 1)
    scale_factor = sdata.images[hires_key].transform[blocco_key].scale[1]
    # Infer geojson path: assumes file is named 'tissue_hires_image_{blocco_key}.geojson' in '~/data/geojson_dir/'
    geojson_path = f"{geojson_dir}tissue_hires_image_{blocco_key}.geojson"
    if not os.path.exists(geojson_path):
        raise FileNotFoundError(f"GeoJSON file not found: {geojson_path}")
    # Read intissue polygons
    intissue_poly = gpd.read_file(geojson_path)
    intissue_scaled = intissue_poly.set_crs(None, allow_override=True)
    # Apply scaling to all geometries in intissue_poly
    intissue_scaled['geometry'] = intissue_scaled['geometry'].apply(
        lambda geom: scale(geom, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))
    )
    
    # Add intissue poly in the sdata
    intissue_rename = "intissue_poly"
    intissue_parse = ShapesModel.parse(intissue_scaled, transformations = {blocco_key: Identity()})
    sdata.shapes[intissue_rename] = intissue_parse
    
    # Extract bins shapes keeping the index 'location_id' for the filtering
    bins_key = [key for key in sdata.shapes if re.search(r'_square_002um$', key)][0]
    sdata[bins_key] = sdata[bins_key].reset_index()  # location_id becomes a column
    # Filter bins intissue
    bins_intissue = sopa.spatial.sjoin(sdata,
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
    sdata.shapes[bins_shape_rename] = bins_intissue # ready to get inside the sdata 


    # Annotate the table with the filtered bins
    sdata["square_002um"].obs["region"] = pd.Categorical([bins_shape_rename] * len(sdata["square_002um"]))
    sdata["square_002um"].uns["spatialdata_attrs"] = {
        "region": bins_shape_rename,  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
    }
    sdata.set_table_annotates_spatialelement("square_002um", region=bins_shape_rename)    
    
    # Filter the table
    sdata['filtered'] = sd.match_table_to_element(sdata, 
        element_name=bins_shape_rename, 
        table_name="square_002um"
    )

    # Map the 'exp_condition' in the 'filtered' table
    location_to_name = sdata[bins_shape_rename].set_index('location_id')['name']
    # Map the names to adata.obs using location_id
    sdata['filtered'].obs['sample_id'] = sdata['filtered'].obs['location_id'].map(location_to_name)

    # create a new sdata with only the interesting elements.
    final = sdata.subset([bins_shape_rename, intissue_rename, 'filtered'], filter_tables = False)
    image_rename = "full_image"
    final[image_rename] = sdata[f'{blocco_key}_full_image'].copy()
    
    return final
  

#-------------------------------------------------------------------------------

def divide_samples(samples_dict, 
                   input_dir  = "/mnt/europa/valerio/data/zarr_store/filtered/",
                   output_dir = "/mnt/europa/valerio/data/zarr_store/samples/"
):
    '''
    The function need a dictionary organized like that:
      block1:
              sample01:
                      min_coordinate = [xmin ,ymin]
                      max_coordinate = [xmax, ymax]
              sample02:
                      min_coordinate = [xmin ,ymin]
                      max_coordinate = [xmax, ymax]
      block2:
              sample03:
                      min_coordinate = [xmin ,ymin]
                      max_coordinate = [xmax, ymax]
    '''
    for blocco, samples in samples_dict.items():
        sdata = read_zarr_standardized(f"{input_dir}filtered_{blocco}.zarr")
        if sdata is None:
            print(f"Warning: {blocco} not found in directory; skipping.")
            continue
        for sample, bbox in samples.items():
            min_coordinate = bbox['min_coordinate']
            max_coordinate = bbox['max_coordinate']
            # Bounding box query
            sdata_bbox = sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coordinate,
                max_coordinate=max_coordinate,
                target_coordinate_system=blocco
            )
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
            
            # renaming coordinate system
            sdata.rename_coordinate_systems(blocco : f"{blocco}_{sample}")

            # Sanity check
            try:
                sanity = sanity_check(sdata_bbox)
                if sanity is None:
                    out_path = f"{output_dir}{blocco}_{sample}.zarr"
                    sdata_bbox.write(out_path)
                    print(f"Wrote {out_path}")
                else:
                    print(f"Sanity check returned a value for {blocco}_{sample}; skipping write.")
            except AssertionError as e:
                print(f"Sanity check failed for {blocco}_{sample}: {e}")
    return print("Job done!")  


