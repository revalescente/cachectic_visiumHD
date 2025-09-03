import re
import os
import sopa
import pandas as pd
import geopandas as gpd
from shapely.affinity import scale
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity
from sopa.io.standardize import sanity_check, write_standardized, read_zarr_standardized

def sdata_pp(path=None):
    # Infer blocco name from path
    # e.g. '/mnt/europa/valerio/data/spe_blocco1.zarr' -> 'blocco1'
    blocco_match = re.search(r'(blocco\d+)', path)
    if not blocco_match:
        raise ValueError("Cannot infer blocco name from path.")
    blocco_key = blocco_match.group(1)

    # Read spatialdata object
    spe = sd.read_zarr(path)

    # Find hires image key using regex
    hires_keys = [key for key in spe.images if re.search(r'_hires_image$', key)]
    if not hires_keys:
        raise ValueError("No hires image found in spe.images")
    hires_key = hires_keys[0]

    # Access the y scale value (index 1)
    scale_factor = spe.images[hires_key].transform[blocco_key].scale[1]
    # Infer geojson path: assumes file is named 'tissue_hires_image_{blocco_key}.geojson' in '~/data/geojson_dir/'
    geojson_path = os.path.expanduser(
        f"~/data/geojson_dir/tissue_hires_image_{blocco_key}.geojson"
    )
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
    bins = spe[bins_shape_name].reset_index()  # location_id becomes a column
    # Filter bins intissue
    bins_intissue = gpd.sjoin(
        bins,
        intissue_scaled[['geometry', 'name']],
        how='inner',
        predicate='intersects'
    )
    # there can be duplicated bins because inside 2 different tissue, I exclude them
    bins_intissue = bins_intissue[bins_intissue['location_id'].duplicated(keep=False) == False]
    # Add in the spe object both the bins and the intissue poly
    bins_shape_rename = f"{blocco_key}_intissue_002um"
    intissue_rename = f"{blocco_key}_intissue"
    intissue_parse = ShapesModel.parse(intissue_scaled, transformations = {blocco_key: Identity()})
    intersection_parse = ShapesModel.parse(bins_intissue, transformations = {blocco_key: Identity()})
    spe.shapes[bins_shape_rename] = intersection_parse
    spe.shapes[intissue_rename] = intissue_parse

    # Annotate the table (just add the region columns to the obs of the table)
    spe["square_002um"].obs["region"] = pd.Categorical([bins_shape_rename] * len(spe["square_002um"]))
    spe["square_002um"].obs[["region", "location_id"]]
    spe["square_002um"].uns["spatialdata_attrs"] = {
        "region": bins_shape_rename,  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
    }

    # Filter the table
    spe['filtered'] = sd.match_table_to_element(spe, element_name=bins_shape_rename, table_name="square_002um")

    # Map the 'exp_condition' in the 'filtered' table
    location_to_name = spe[bins_shape_rename].set_index('location_id')['name']
    # Map the names to adata.obs using location_id
    spe['filtered'].obs['exp_cond'] = spe['filtered'].obs['location_id'].map(location_to_name)

    # New column for sample identification blocco_expcond
    spe['filtered'].obs['sample_id'] = f"{blocco_key}_" + spe['filtered'].obs['exp_cond'].astype(str)
    
    # create a new sdata with only the interesting elements.
    full_image_name = f"{blocco_key}_full_image"
    final = spe.subset([full_image_name, bins_shape_rename, intissue_rename, 'filtered'], filter_tables = False)
    return final

#-------------------------------------------------------------------------------

def divide_samples(spe_blocks, blocco_sample_bbox_dict, output_dir="/mnt/europa/valerio/data/zarr_store/blocchi/"):
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
    for blocco, samples in blocco_sample_bbox_dict.items():
        spe = spe_blocks.get(blocco)
        if spe is None:
            print(f"Warning: {blocco} not found in spe_blocks; skipping.")
            continue
        for sample, bbox in samples.items():
            min_coordinate = bbox['min_coordinate']
            max_coordinate = bbox['max_coordinate']
            # Bounding box query
            sdata_bbox = spe.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coordinate,
                max_coordinate=max_coordinate,
                target_coordinate_system=blocco,
            )
            # Subset elements
            subset_keys = [f"{blocco}_full_image", f"{blocco}_intissue_002um", f"{blocco}_intissue", "filtered"]
            sdata_subset = sdata_bbox.subset(subset_keys, filter_tables=False)
            
            # --- INTEGRATION: filter bins and table ---
            # 1. Filter 'blocco_intissue_002um' by exp_condition (column "name")
            element_key = f"{blocco}_intissue_002um"
            table_key = "filtered"
            if element_key in sdata_subset:
                sdata_subset[element_key] = sdata_subset[element_key][sdata_subset[element_key]['name'] == sample]
            else:
                print(f"Element {element_key} not found in sdata_subset for {blocco}_{sample}")

            # 2. Update 'filtered' table to match filtered element
            if hasattr(sd, "match_table_to_element"):  # sd must be your SpatialData module
                if table_key in sdata_subset:
                    sdata_subset[table_key] = sd.match_table_to_element(sdata_subset, element_name=element_key, table_name=table_key)
                else:
                    print(f"Table {table_key} not found in sdata_subset for {blocco}_{sample}")
            else:
                print("sd.match_table_to_element not found. Make sure sd (SpatialData) is imported.")

            # sopa metadata 
            sopa.utils.set_sopa_attrs(
                sdata_subset,
                cell_segmentation_key = f"{blocco}_full_image",
                tissue_segmentation_key = f"{blocco}_full_image",
                bins_table_key = "filtered",
                boundaries_key = f"{blocco}_intissue"
            )
            
            # Sanity check
            try:
                sanity = sanity_check(sdata_subset)
                if sanity is None:
                    out_path = f"{output_dir}{blocco}_{sample}.zarr"
                    write_standardized(sdata_subset, out_path)
                    print(f"Wrote {out_path}")
                else:
                    print(f"Sanity check returned a value for {blocco}_{sample}; skipping write.")
            except AssertionError as e:
                print(f"Sanity check failed for {blocco}_{sample}: {e}")


