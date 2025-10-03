import re
import os
import sopa
import pandas as pd
import geopandas as gpd
from shapely.affinity import scale
import spatialdata as sd
from spatialdata.transformations import Identity
from sopa.io.standardize import sanity_check, read_zarr_standardized
from spatialdata_io import visium_hd
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from spatialdata.transformations import Identity
from spatialdata import SpatialData
import json

def preprocess_step(block_numbers = None, 
                       input_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out/",
                       output_path = "/mnt/europa/valerio/data/zarr_store/",
                       images_path = "/mnt/europa/valerio/HE_images/color_corrected/"
):
    """
    Function to prepare data to segmentation: filtering useless things
    
    if you don't specify which blocks to process the function will process all of them
    
    """
    # Which block to preprocess
    if block_numbers is None:
      block_numbers = [1,2,3,4,5,6,7,9]
    else:
      block_numbers

    sdata_blocks = {}
    
    # 1. reading data from source 
    for i in block_numbers:
        blocco_key = f"blocco{i}"
        sdata_blocks[block_name] = visium_hd(
            path=f"{input_path}blocco{i}/outs",
            dataset_id=f"blocco{i}",
            filtered_counts_file=False,
            bin_size='002',
            bins_as_squares=True,
            annotate_table_by_labels=False,
            fullres_image_file=f"{images_path}pp_blocco{i}_20x.tif",
            load_all_images=False,
            var_names_make_unique=True
        )
    # 1a. Writing data as zarr stores.
    for blocco_key, sdata in sdata_blocks.items():
      sdata.write(f"{output_path}general/{block_key}.zarr")
    
    # 2. Filtering of every blocks to remove useless bins and other little things
    for num in block_numbers:
      path = f'/mnt/europa/valerio/data/zarr_store/general/blocco{num}.zarr'
      try:
          result = pp.sdata_pp(path)
          print(f"Processed blocco {num}")
          # resaving as filtered data
          result.write(f'{output_path}filtered/filtered_blocco{num}.zarr')
          print(f"Zarr Store created for blocco{num}")
      except Exception as e:
          print(f"Error processing blocco {num}: {e}")
      return print("Job done")

# ------------------------------------------------------------------------------

def sdata_pp(path = None
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

    # Extract bins shapes keeping the index 'location_id' for the filtering
    bins_shape_name = f"{blocco_key}_square_002um"
    bins = sdata[bins_shape_name].reset_index()  # location_id becomes a column
    # Filter bins intissue
    bins_intissue = gpd.sjoin(
        bins,
        intissue_scaled[['geometry', 'name']],
        how='inner',
        predicate='intersects'
    )
    # there can be duplicated bins because inside 2 different tissue, I exclude them
    bins_intissue = bins_intissue[bins_intissue['location_id'].duplicated(keep=False) == False]
    # Add in the sdata object both the bins and the intissue poly
    bins_shape_rename = "intissue_002um"
    intissue_rename = "intissue_poly"
    intissue_parse = ShapesModel.parse(intissue_scaled, transformations = {blocco_key: Identity()})
    intersection_parse = ShapesModel.parse(bins_intissue, transformations = {blocco_key: Identity()})
    sdata.shapes[bins_shape_rename] = intersection_parse
    sdata.shapes[intissue_rename] = intissue_parse

    # Annotate the table (just add the region columns to the obs of the table)
    sdata["square_002um"].obs["region"] = pd.Categorical([bins_shape_rename] * len(sdata["square_002um"]))
    sdata["square_002um"].uns["spatialdata_attrs"] = {
        "region": bins_shape_rename,  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
    }

    # Filter the table
    sdata['filtered'] = sd.match_table_to_element(sdata, element_name=bins_shape_rename, table_name="square_002um")

    # Map the 'exp_condition' in the 'filtered' table
    location_to_name = sdata[bins_shape_rename].set_index('location_id')['name']
    # Map the names to adata.obs using location_id
    sdata['filtered'].obs['exp_cond'] = sdata['filtered'].obs['location_id'].map(location_to_name)

    # create a new sdata with only the interesting elements.
    full_image_name = "full_image"
    final = sdata.subset([full_image_name, bins_shape_rename, intissue_rename, 'filtered'], filter_tables = False)
    return final

#-------------------------------------------------------------------------------

def divide_samples(samples_dict, 
                   input_dir  = "/mnt/europa/valerio/data/zarr_store/filtered/",
                   output_dir = "/mnt/europa/valerio/data/zarr_store/blocchi/"
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
            # Subset elements
            subset_keys = ["full_image", "intissue_002um", "intissue", "filtered"]
            sdata_subset = sdata_bbox.subset(subset_keys, filter_tables=False)
            
            # --- INTEGRATION: filter bins, intissue and table ---
            # 1. Filter 'blocco_intissue_002um' and 'blocco_intissue' by exp_condition (column "name") 
            table_key = "filtered"
            
            # 1a. filtering bins not belonging to the right sample
            bins_key = "intissue_002um"
            if bins_key in sdata_subset:
                sdata_subset[bins_key] = sdata_subset[bins_key][sdata_subset[bins_key]['name'] == sample]
            else:
                print(f"Element {bins_key} not found in sdata_subset for {blocco}_{sample}")
            
            # 1b. filtering intissue poly for the same reason
            intissue_key = "intissue"
            if intissue_key in sdata_subset:
                sdata_subset[intissue_key] = sdata_subset[intissue_key][sdata_subset[intissue_key]['name'] == sample]
            else:
                print(f"Element {intissue_key} not found in sdata_subset for {blocco}_{sample}")

            # 2. Update 'filtered' table to match filtered bins 
            if hasattr(sd, "match_table_to_element"):  # sd must be your SpatialData module
                if table_key in sdata_subset:
                    sdata_subset[table_key] = sd.match_table_to_element(sdata_subset, element_name=bins_key, table_name=table_key)
                else:
                    print(f"Table {table_key} not found in sdata_subset for {blocco}_{sample}")
            else:
                print("sd.match_table_to_element not found. Make sure sd (SpatialData) is imported.")

            # sopa metadata 
            sopa.utils.set_sopa_attrs(
                sdata_subset,
                cell_segmentation_key = "full_image",
                tissue_segmentation_key = "full_image",
                bins_table_key = "filtered"
            )

            # Sanity check
            try:
                sanity = sanity_check(sdata_subset)
                if sanity is None:
                    out_path = f"{output_dir}{blocco}_{sample}.zarr"
                    sdata_subset.write(out_path)
                    print(f"Wrote {out_path}")
                else:
                    print(f"Sanity check returned a value for {blocco}_{sample}; skipping write.")
            except AssertionError as e:
                print(f"Sanity check failed for {blocco}_{sample}: {e}")


