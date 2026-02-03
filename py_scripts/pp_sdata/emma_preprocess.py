import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import re
import os
import geopandas as gpd
import pandas as pd
from shapely.affinity import scale
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp

geojson_dir = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
blocco_key = 'blocco2'
sample_key = 'c26'
hires_key = f'{blocco_key}_hires_image'
geojson_path = f"{geojson_dir}/{blocco_key}_{sample_key}.geojson"

sdata = visium_hd(
        path= f"/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/{blocco_key}/outs",
        dataset_id=f"{blocco_key}",
        filtered_counts_file=False,
        bin_size=['008','016'],
        bins_as_squares=True,
        annotate_table_by_labels=False,
        load_all_images=False,
        var_names_make_unique=True,
        image_models_kwargs = {'dims' : ['c', 'y', 'x']}
    )

# # Find hires image key using regex
# hires_keys = [key for key in sdata.images if re.search(r'_hires_image$', key)]
# hires_key = hires_keys[0]
# 
# # Access the y scale value (index 1)
# scale_factor = sdata.images[hires_key].transform[blocco_key].scale[1]
# # Read intissue polygons
# intissue_poly = gpd.read_file(geojson_path)
# intissue_scaled = intissue_poly.set_crs(None, allow_override=True)
# # Apply scaling to all geometries in intissue_poly
# intissue_scaled['geometry'] = intissue_scaled['geometry'].apply(
#     lambda geom: scale(geom, xfact=scale_factor, yfact=scale_factor, origin=(0, 0))
# )

#blocco 2 intissue diverso
geojson_path = f"{geojson_dir}in_tissue/tissue_fullres_image_{blocco_key}.geojson"
intissue_poly = gpd.read_file(geojson_path)
intissue_poly = intissue_poly.set_crs(None, allow_override=True)
# Add intissue poly in the sdata
intissue_rename = "intissue_poly"
intissue_parse = ShapesModel.parse(intissue_poly, transformations = {blocco_key: Identity()})
sdata.shapes[intissue_rename] = intissue_parse
    
# Extract bins shapes keeping the index 'location_id' for the filtering
bins_key = [key for key in sdata.shapes if re.search(r'_square_008um$', key)][0]
sdata[bins_key] = sdata[bins_key].reset_index()  # location_id becomes a column
# Filter bins intissue
bins_intissue = sopa.spatial.sjoin(sdata,
      bins_key,
      intissue_rename,
      how = "inner",
      predicate = "intersects",
      target_coordinate_system = blocco_key
)

# remove duplicated
bins_intissue = bins_intissue[bins_intissue['location_id'].duplicated(keep=False) == False]
bins_intissue = bins_intissue[['location_id', 'geometry', 'name']]
# Add in the sdata object both the bins and the intissue poly
bins_shape_rename = "intissue_008um"
sdata.shapes[bins_shape_rename] = bins_intissue # ready to get inside the sdata 

# Annotate the table with the filtered bins
sdata["square_008um"].obs["region"] = pd.Categorical([bins_shape_rename] * len(sdata["square_008um"]))
sdata["square_008um"].uns["spatialdata_attrs"] = {
        "region": bins_shape_rename,  # name of the Shapes element we will use later
        "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
        "instance_key": "location_id",  # column that matches a given obs in the table to a given circle
}
sdata.set_table_annotates_spatialelement("square_008um", region=bins_shape_rename)    
    
# Filter the table
sdata['filtered'] = sd.match_table_to_element(sdata, 
        element_name=bins_shape_rename, 
        table_name="square_008um"
)

# Map the 'exp_condition' in the 'filtered' table
location_to_name = sdata[bins_shape_rename].set_index('location_id')['name']
# Map the names to adata.obs using location_id
sdata['filtered'].obs['sample_id'] = sdata['filtered'].obs['location_id'].map(location_to_name)

# Map the 'exp_condition' in the 'filtered' table
location_to_name = sdata[bins_shape_rename].set_index('location_id')['name']
# Map the names to adata.obs using location_id
sdata['filtered'].obs['sample_id'] = sdata['filtered'].obs['location_id'].map(location_to_name)

# GFP polys --------------------------------------------------------------------
gfp_path = f"/mnt/europa/valerio/data/json/geojson_dir/GFP_inflamed/{blocco_key}_GFP_and_inflamed.geojson"
gfp_poly = gpd.read_file(gfp_path)
gfp_poly = gfp_poly.set_crs(None, allow_override=True)
# proviamo a separare i polygoni cosi matplot fa meno fatica... ???
#gfp_poly = gfp_poly.explode(index_parts=False).reset_index(drop=True)
gfp_parse = ShapesModel.parse(gfp_poly, transformations={blocco_key: Identity()})

sdata.shapes["GFP_poly"] = gfp_parse

table_name = 'filtered'
table = sdata[table_name].copy()
attrs = table.uns['spatialdata_attrs']
region_name = attrs['region']
instance_key = attrs['instance_key']
# Initialize annotation column
table.obs['in_treatment'] = False
# 2. Determine Coordinate System (pick the first available)
bins_element = sdata[region_name].copy()
coord_system = list(get_transformation(bins_element, get_all=True).keys())[0]
# 3. Perform SINGLE Spatial Join
# Identifies all bins overlapping ANY polygon in 'GFP_poly'
joined_bins = sopa.spatial.sjoin(
    sdata, 
    region_name, 
    'GFP_poly', 
    how='inner', 
    predicate='intersects', 
    target_coordinate_system=coord_system
)
if not joined_bins.empty:
    # 4. Identify IDs
    # Group A: IDs to REMOVE (infiammazione)
    mask_inf = joined_bins['name_right'].str.contains("infiammazione", case=False, na=False)
    ids_to_remove = joined_bins[mask_inf].index.unique()
    # Group B: IDs to ANNOTATE (fibre_trattate)
    mask_treat = joined_bins['name_right'].str.contains("fibre_trattate", case=False, na=False)
    ids_to_treat = joined_bins[mask_treat].index.unique()
    # Prevent annotating bins that are about to be removed
    ids_to_treat = ids_to_treat.difference(ids_to_remove)
    # 5. Update the Table
    # First, apply annotations to the table (before filtering rows)
    if len(ids_to_treat) > 0:
        is_treated = table.obs[instance_key].isin(ids_to_treat)
        table.obs.loc[is_treated, 'in_treatment'] = True
    # Second, filter the table rows (remove inflammation)
    if len(ids_to_remove) > 0:
        keep_mask = ~table.obs[instance_key].isin(ids_to_remove)
        table = table[keep_mask].copy()
    sdata[table_name] = table
    # 6. Sync Shapes to Table
    # match_element_to_table returns a tuple: ({element_name: element}, table)
    filtered_elements, _ = sd.match_element_to_table(sdata, region_name, table_name)
    # Extract the specific element (GeoDataFrame) from the dictionary and assign it
    sdata[region_name] = filtered_elements[region_name]
    print(f"Processed {region_name}: Removed {len(ids_to_remove)} bins. Annotated {len(ids_to_treat)} bins.")
else:
    print(f"Processed {region_name}: No overlap with polygons.")

# grafico di controllo ---------------------------------------------------------
plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_images(f'{blocco_key}_hires_image', scale = 'scale3'
).pl.render_shapes('intissue_008um', outline=False, outline_alpha=1, outline_width=1, fill_alpha=1
).pl.show(ax = ax, coordinate_systems=blocco_key, 
    save = f'output_python/roba_da_buttare/{blocco_key}_emma.png'
)

# separazione campioni ---------------------------------------------------------

samples_dict = pp.read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')
block_samples = samples_dict.get(blocco_key)
output_dir = "/mnt/europa/shared/sandri_visiumHD_data/bins/"


if not block_samples:
    print(f"No samples found for {blocco_key} in the dictionary.")
else:
    print(f"Processing samples for {blocco_key}...")
    for sample_name, sample_info in block_samples.items():
        print(f"\n--- Processing sample: {sample_name} ---")
        min_coord = sample_info['min_coordinate']
        max_coord = sample_info['max_coordinate']
        print(f"  BBox: Min {min_coord}, Max {max_coord}")
        # 3. Bounding Box Query
        try:
            sdata_bbox = sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=blocco_key
            )
        except Exception as e:
            print(f"  Error during bounding box query: {e}")
            continue
        # 4. Filter the 'filtered' table by sample_id
        # Note: You used 'sample_id' in get_values, ensuring this column exists in table.obs
        try:
            # Check if table exists in the cropped object
            if "filtered" not in sdata_bbox.tables:
                print("  'filtered' table not found in bbox crop. Skipping.")
                continue
            table = sdata_bbox["filtered"]
            target_id = sample_name # e.g., 'c26STAT3'
            subset_mask = table.obs['sample_id'] == target_id
            sdata_bbox["filtered"] = table[subset_mask].copy()
            print(f"  Filtered table to {target_id}: {sdata_bbox['filtered'].n_obs} bins")
        except KeyError:
            print(f"  'sample_id' column not found in table.obs. Skipping table filtering.")
        except Exception as e:
            print(f"  Error filtering table: {e}")
            continue
        # 5. Match Elements (Shapes) to the Filtered Table
        # Your example used 'intissue_002um'. I'll generalize to what's in your sdata (intissue_008um)
        # Adjust 'element_name' to matches what is actually in your sdata structure
        target_shape_key = 'intissue_008um' # Based on your sdata print output
        if target_shape_key in sdata_bbox.shapes:
            try:
                matched_elements = sd.match_element_to_table(
                    sdata_bbox, 
                    element_name=target_shape_key, 
                    table_name='filtered'
                )
                # match_element_to_table returns (elements_dict, table)
                sdata_bbox[target_shape_key] = matched_elements[0][target_shape_key]
                print(f"  Matched {target_shape_key} to table.")
            except Exception as e:
                 print(f"  Error matching {target_shape_key}: {e}")
        # 6. Filter 'intissue_poly' (Annotated Regions)
        if 'intissue_poly' in sdata_bbox.shapes:
            try:
                poly_gdf = sdata_bbox['intissue_poly']
                # Filtering by 'name' column matching the sample
                sdata_bbox['intissue_poly'] = poly_gdf[poly_gdf['name'] == sample_name]
                print(f"  Filtered 'intissue_poly' for {sample_name}")
            except Exception as e:
                print(f"  Error filtering intissue_poly: {e}")
        # 7. Setup Sopa Metadata
        try:
            # Check if keys exist before setting attributes to avoid errors
            # Your manual code used "full_image" which might not be in the crop if the crop 
            # logic renamed it or if it wasn't selected. 
            # Bounding box query usually keeps image names.
            # Assuming 'blocco1_hires_image' is the one you want, based on your sdata print
            # Or pass the generic key if you have a standardized pipeline
            img_key = f"{blocco_key}_hires_image" 
            sopa.utils.set_sopa_attrs(
                sdata_bbox,
                cell_segmentation_key = img_key,   # Updated to match your sdata
                tissue_segmentation_key = img_key, # Updated to match your sdata
                bins_table_key = "filtered"
            )
        except Exception as e:
             print(f"  Warning: Could not set sopa attributes: {e}")
        # 8. Write to Disk
        out_filename = f"{blocco_key}_{sample_name}"
        out_path = os.path.join(output_dir, out_filename)
        try:
            # Optional: Add your sanity_check() here if defined
            # sanity = sanity_check(sdata_bbox)
            print(f"  Writing to {out_path}...")
            sdata_bbox.write(out_path)
            print(f"  Successfully wrote {out_filename}")
        except Exception as e:
            print(f"  Failed to write {out_filename}: {e}")
