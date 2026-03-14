import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import re
import os
import geopandas as gpd
import pandas as pd
from skimage import io
from shapely.affinity import scale
from spatialdata.models import ShapesModel, Image2DModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp

# keys name
BLOCCO_KEY = 'blocco3'
SAMPLE_KEY = 'c26murf1'
HIRES_KEY = f'{BLOCCO_KEY}_hires_image'
FULL_IMAGE_KEY = 'full_image'
FLUO_KEY = 'fluo_image'

# paths of interest
intissue_gfp_dir = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
data_blocco2 = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell"
data_b_all = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_cellexpans"
arivis_dir = "/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks"
fullres_path = "/mnt/europa/valerio/HE_images/color_corrected/samples"
fluo_path = "/mnt/europa/valerio/Fluo_images/warped_tif/samples"



# 1. load the data
sdata = visium_hd(
        path= f"{data_b_all}/{BLOCCO_KEY}/outs",
        dataset_id=f"{BLOCCO_KEY}",
        filtered_counts_file=False,
        bin_size=['002','008','016'],
        bins_as_squares=True,
        annotate_table_by_labels=False,
        load_all_images=False,
        var_names_make_unique=True,
        image_models_kwargs = {'dims' : ['c', 'y', 'x']},
        load_segmentations_only = False,
        load_nucleus_segmentations = False
)

# 1b. remove the cell segmentations of spaceranger
del sdata[f'{BLOCCO_KEY}_cell_segmentations']
del sdata['cell_segmentations']

# 2. separate the 3 tissues
samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')
block_samples = samples_dict.get(BLOCCO_KEY)

if not block_samples:
    print(f"No samples found for {BLOCCO_KEY} in the dictionary.")
else:
    print(f"Processing samples for {BLOCCO_KEY}...")
    # Iterate over each sample of the blocco
    for sample_name, sample_info in block_samples.items():
        print(f"\n--- Processing sample: {sample_name} ---")
        min_coord = sample_info['min_coordinate']
        max_coord = sample_info['max_coordinate']
        print(f"  BBox: Min {min_coord}, Max {max_coord}")
        
        # Bounding Box Query
        try:
            sdata_bbox = sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=BLOCCO_KEY
            )
        except Exception as e:
            print(f"  Error during bounding box query: {e}")
            continue
        
        # 2. ADD FULL RES IMAGE
        full_image = io.read_image
        
        # 2a. Add intissue_GPF_polys to the sdata_bbox
        all_poly = gpd.read_file(intissue_gfp_dir + f"/{BLOCCO_KEY}_{sample_name}.geojson")
        all_poly = all_poly.set_crs(None, allow_override=True)
        
        # Store the full parsed polys inside the sdata_bbox shapes
        polys_name = f"intissue_GFP_poly_{sample_name}"
        polys_parse = ShapesModel.parse(all_poly, transformations={BLOCCO_KEY: Identity()})
        sdata_bbox.shapes[polys_name] = polys_parse

        # 2b. Spatial Join: Annotate the bins/tables with True/False for each geometry
        bin_keys = [f'{BLOCCO_KEY}_square_002um', f'{BLOCCO_KEY}_square_008um', f'{BLOCCO_KEY}_square_016um']
        
        for idx, row in all_poly.iterrows():
            geom_name = row['name']  # e.g., 'c26SMAD23', 'infiammazione_gr1', 'fibre_trattate_gr1'
            
            # Create a GeoDataFrame for this single geometry
            single_poly_gdf = gpd.GeoDataFrame([row], geometry='geometry', crs=all_poly.crs)
            
            # Iterate through the bin sizes and apply spatial join
            for bin_key in bin_keys:
                if bin_key in sdata_bbox.shapes:
                    table_name = bin_key.replace(f'{BLOCCO_KEY}_', '') # e.g., 'square_002um'
                    
                    if table_name in sdata_bbox.tables:
                        # Fetch the bins shapes for the current resolution
                        bins_gdf = sdata_bbox.shapes[bin_key]
                        
                        # Use geopandas spatial join to find intersections
                        # (this mimics sopa.spatial.sjoin but gives you direct True/False boolean arrays)
                        joined = gpd.sjoin(bins_gdf, single_poly_gdf, how='left', predicate='intersects')
                        
                        # 'index_right' will not be NA for bins that intersect with the geometry
                        intersecting_indices = joined[joined['index_right'].notna()].index
                        
                        # Assign boolean True/False back to the table's .obs
                        # Initialize column as False
                        obs_col_name = f"in_{geom_name}"
                        sdata_bbox.tables[table_name].obs[obs_col_name] = False
                        
                        # Update intersecting rows to True
                        # Ensure we only update indices that actually exist in the table's .obs
                        valid_intersecting = intersecting_indices.intersection(sdata_bbox.tables[table_name].obs.index)
                        sdata_bbox.tables[table_name].obs.loc[valid_intersecting, obs_col_name] = True
                        
                        print(f"    -> Added boolean mask '{obs_col_name}' to table '{table_name}' ({len(valid_intersecting)} True)")
