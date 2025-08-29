import re
import os
import pandas as pd
import geopandas as gpd
from shapely.affinity import scale
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity

def sdata_pp(path='/mnt/europa/valerio/data/spe_blocco1.zarr'):
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

    # Add in the spe object
    intersection_parse = ShapesModel.parse(bins_intissue, transformations={blocco_key: Identity()})
    spe.shapes['intissue_002um'] = intersection_parse

    # Annotate the table (just add the region columns to the obs of the table)
    table_name = "square_002um"
    spe[table_name].obs["region"] = pd.Categorical(["intissue_002um"] * len(spe[table_name]))
    spe[table_name].obs[["region", "location_id"]]

    spe[table_name].uns["spatialdata_attrs"] = {
        "region": "intissue_002um",  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
    }

    # Filter the table
    spe['filtered'] = sd.match_table_to_element(spe, element_name='intissue_002um', table_name=table_name)

    full_image_name = f"{blocco_key}_full_image"
    final = spe.subset([full_image_name, 'intissue_002um', 'filtered'], filter_tables = False)
    return final
