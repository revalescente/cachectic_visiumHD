import spatialdata as spd
import gc
import geopandas as gpd
import py_scripts.utils.spaceranger_utility as sr
from py_scripts.utils.utils_fun import read_from_json
from spatialdata.models import ShapesModel, Image2DModel
from spatialdata.transformations import Identity, set_transformation, get_transformation, Scale
from skimage import io
import sopa
import numpy as np

# spaceranger pipeline update
blocco = "blocco1"
# read sdata full blocco
sdata = spd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/blocchi/{blocco}")

save_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/"
poly_path = "/mnt/europa/valerio/data/json/geojson_dir/in_tissue/fullres/"
fullres_path = "/mnt/europa/valerio/HE_images/color_corrected/"
fluo_path = "/mnt/europa/valerio/Fluo_images/warped_tif/"
json_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell/"

"""
Iteratively creates and saves cropped subsets of a SpatialData object,
following the user's specified procedural script logic for each tissue type.
Args:
    sdata (spd.SpatialData): The input SpatialData object.
    save_path (str): The base directory where the output .zarr files will be saved.
"""
# --- 1. Initial Setup outside the loop ---
blocco_key = sdata.path.stem
table_key = "segmentation_counts"
coord_key = "downscale_to_hires"
fullres_coord_key = sdata.path.stem

# Define original element keys
intissue_key_orig = f"{blocco_key}_intissue_boundaries"
nuclei_key_orig = f"{blocco_key}_nuclei_boundaries"
all_poly_key_orig = f"{blocco_key}_all_poly"
fullimage_key_orig = f"{blocco_key}_full_image"
fluoimage_key_orig = f"{blocco_key}_GFP_image"

# Add the 'all_poly
poly = gpd.read_file(f"{poly_path}{blocco_key}_allpolys_fullres.geojson")
poly = poly.set_crs(None, allow_override=True)
poly_parse = ShapesModel.parse(poly, transformations={fullres_coord_key: Identity()})
sdata[all_poly_key_orig] = poly_parse

# Add the fullres image
img = io.imread(f"{fullres_path}pp_{blocco_key}_20x.tif")
img_parsed = Image2DModel.parse(data=img, 
            scale_factors=(2, 2, 2), 
            transformations={fullres_coord_key: Identity()},
            dims=("y", "x", "c")
) 
sdata[fullimage_key_orig] = img_parsed

# Add the fluo image
fluo = io.imread(f"{fluo_path}{blocco_key}.tif")

fluo_rescaled = np.empty_like(fluo)
for c in range(fluo.shape[-1]):
    ch = fluo[..., c]
    ch_min = ch.min()
    ch_max = ch.max()
    # Avoid division by zero if channel is flat
    if ch_max > ch_min:
        fluo_rescaled[..., c] = 255 * (ch - ch_min) / (ch_max - ch_min)
    else:
        fluo_rescaled[..., c] = 0  # or ch_min, or 255, as appropriate

fluo_parsed = Image2DModel.parse(data=fluo_rescaled, 
            scale_factors=(2, 2, 2), 
            transformations={fullres_coord_key: Identity()},
            dims=("y", "x", "c")
) 
sdata[fluoimage_key_orig] = fluo_parsed

# Add transformation to the nuclei boundaries
scalejs = read_from_json(f'{json_path}{blocco_key}/outs/segmented_outputs/spatial/scalefactors_json.json')
scale_factor = 1 / scalejs['tissue_hires_scalef']
scaler_transform = Scale([scale_factor, scale_factor], axes=("x", "y"))
set_transformation(
    element=sdata[nuclei_key_orig],
    transformation=scaler_transform,
    to_coordinate_system=fullres_coord_key
)
set_transformation(
    element=sdata[intissue_key_orig],
    transformation=scaler_transform,
    to_coordinate_system=fullres_coord_key
)

# Add GFP channel info to the table
channel_aggregation = sopa.aggregation.aggregate_channels(
  sdata, image_key=fluoimage_key_orig, shapes_key=nuclei_key_orig, 
  expand_radius_ratio=0, mode='max', no_overlap=False
) 
max_values_vector = channel_aggregation.max(axis=1)
sdata[table_key].obs['GFP_value'] = max_values_vector

# old blocco 2 problem, maybe solved
# Get all unique tissue types to iterate over
#sdata.tables[table_key].obs['tissue_type'] = sdata.tables[table_key].obs['tissue_type'].replace('c26foxO', 'c26murf1')
#
# # 2) If it's a categorical column, you need to also update the categories
# if sdata.tables[table_key].obs['tissue_type'].dtype.name == 'category':
#     sdata.tables[table_key].obs['tissue_type'] = sdata.tables[table_key].obs['tissue_type'].cat.rename_categories({'c26foxO': 'c26murf1'})
#
# # 3) Verify
print(sdata.tables[table_key].obs['tissue_type'].unique())

sdata. shapes["blocco2_all_poly"]["name"] = (
    sdata.shapes["blocco2_all_poly"]["name"]
    .replace('c26foxO', 'c26murf1')
)

# Verify
print(sdata.shapes["blocco2_all_poly"]["name"]. unique())

# sdata.delete_element_from_disk(table_key)
# sdata.write_element(table_key)

all_tissue_types = sdata[table_key].obs['tissue_type'].unique().tolist()
print(f"Found tissue types in '{blocco_key}': {all_tissue_types}")
# --- 2. Loop through each tissue type and apply the script's logic ---
for tissue_type in all_tissue_types:
    print(f"\n--- Processing tissue: {tissue_type} ---")
    # Define keys for the new, filtered elements for this iteration
    nuclei_key_filtered = f"{blocco_key}_{tissue_type}_nuclei_boundaries"
    intissue_key_filtered = f"{blocco_key}_{tissue_type}_intissue_boundaries"
    fullimage_key_filtered = f"{blocco_key}_{tissue_type}_full_image"
    fluoimage_key_filtered = f"{blocco_key}_{tissue_type}_GFP_image"
    gfp_key = "in_treatment_poly"
    inflame_key = "inflamed_poly"
    # a. Filter shapes and add them as NEW elements to the main sdata object
    sdata[intissue_key_filtered] = sdata[all_poly_key_orig][sdata[all_poly_key_orig]['name'] == tissue_type].copy()
    sdata[nuclei_key_filtered] = sdata[nuclei_key_orig][sdata[nuclei_key_orig]['name'] == tissue_type].copy()
    sdata[gfp_key] = sdata[all_poly_key_orig][sdata[all_poly_key_orig]['name'].str.contains('fibre_trattate')]
    sdata[inflame_key] = sdata[all_poly_key_orig][sdata[all_poly_key_orig]['name'].str.contains('infiammazione')]
    # b. Get the extent of the newly created filtered tissue boundary
    extent = spd.get_extent(
        sdata, 
        coordinate_system=fullres_coord_key, 
        exact=True, 
        elements=[intissue_key_filtered]
    )
    # c. Crop the entire sdata object (including the newly added elements)
    cropped_sdata = sr.crop0(sdata, fullres_coord_key, extent)
    # d. Assemble a temporary sdata object from the cropped elements
    
    # Required shapes
    shapes_dict = {
        nuclei_key_filtered:  cropped_sdata[nuclei_key_filtered],
        intissue_key_filtered: cropped_sdata[intissue_key_filtered]
    }
    # Optional shapes
    for key in [gfp_key, inflame_key]:
        if key in cropped_sdata.shapes:
            shapes_dict[key] = cropped_sdata[key]
        else:
            print(f"'{key}' not found, skipping.")
    
    # Assemble tmp_sdata
    # Note: The table here is cropped but not yet filtered by tissue_type
    tmp_sdata = spd.SpatialData(
        images={
            fullimage_key_filtered: cropped_sdata[fullimage_key_orig],
            fluoimage_key_filtered: cropped_sdata[fluoimage_key_orig]
        },
        shapes=shapes_dict,
        tables={table_key: cropped_sdata[table_key]}
    )
    # e. Update the table metadata
    # At this point, the table in tmp_sdata may contain cells from other tissue types
    # that happened to fall within the crop extent. We now filter it.
    cells_to_keep = tmp_sdata[table_key].obs['tissue_type'] == tissue_type
    tmp_sdata.tables[table_key] = tmp_sdata.tables[table_key][cells_to_keep, :].copy()
    # Update region column and link table to the correct shapes
    tmp_sdata[table_key].obs["region"] = nuclei_key_filtered
    tmp_sdata.set_table_annotates_spatialelement(
        table_key, 
        region=nuclei_key_filtered,
        region_key="region",
        instance_key="cell_id"
    )
    # f. Add centroids and sample_id
    centers = spd.get_centroids(tmp_sdata[nuclei_key_filtered], coordinate_system=fullres_coord_key).compute()
    tmp_sdata[table_key].obs['x'] = centers['x']
    tmp_sdata[table_key].obs['y'] = centers['y']
    tmp_sdata[table_key].obs['sample_id'] = f"{blocco_key}_{tissue_type}"
    # ho tutto in linea, ora devo filtrare intissue e poi inserire le colonne .obs per il filtraggio degli altri 2 poligoni 
    if 'index_right' in tmp_sdata[nuclei_key_filtered].columns:
        tmp_sdata[nuclei_key_filtered] = tmp_sdata[nuclei_key_filtered].drop(columns=['index_right'])
        print("Removed 'index_right' column")
    else:
        print("No 'index_right' column found")
    nuclei_intissue = sopa.spatial.sjoin(tmp_sdata,
          nuclei_key_filtered,
          intissue_key_filtered,
          how = "inner",
          predicate = "intersects",
          target_coordinate_system = fullres_coord_key
    )
    nuclei_intissue = nuclei_intissue[['cell_id', 'geometry', 'name_left']]
    nuclei_intissue = nuclei_intissue.rename(columns={'name_left': 'name'})
    # remove duplicated
    nuclei_intissue = nuclei_intissue[nuclei_intissue['cell_id'].duplicated(keep=False) == False]
    # add to tmp_sdata
    tmp_sdata[nuclei_key_filtered] = nuclei_intissue # ready to get inside the sdata 
    # Filter the table
    tmp_sdata[table_key] = spd.match_table_to_element(tmp_sdata, 
            element_name=nuclei_key_filtered, 
            table_name=table_key
    )
    # ------------------------------------------------
    # add info about in_treatment and inflamed
    # ------------------------------------------------
    # --- 1) Initialize annotation columns (always) ---
    tmp_sdata.tables[table_key].obs['in_treatment'] = False
    tmp_sdata.tables[table_key].obs['to_discard'] = False
    # --- 2) Spatial join for IN_TREATMENT (only if gfp_key exists) ---
    if gfp_key in tmp_sdata.shapes:
        joined_treatment = sopa.spatial.sjoin(
            tmp_sdata,
            nuclei_key_filtered,
            gfp_key,
            how='inner',
            predicate='intersects',
            target_coordinate_system=fullres_coord_key
        )
        if not joined_treatment.empty:
            ids_in_treatment = joined_treatment.index.unique()
            mask_treatment = tmp_sdata.tables[table_key].obs.index.isin(ids_in_treatment)
            tmp_sdata.tables[table_key].obs.loc[mask_treatment, 'in_treatment'] = True
            print(f"Marked {mask_treatment.sum()} nuclei as 'in_treatment'")
        else:
            print("No nuclei found inside 'in_treatment_poly'")
    else:
        print(f"Skipping 'in_treatment':  '{gfp_key}' not found in shapes")
    # --- 3) Spatial join for TO_DISCARD (only if inflame_key exists) ---
    if inflame_key in tmp_sdata.shapes:
        joined_inflamed = sopa.spatial.sjoin(
            tmp_sdata,
            nuclei_key_filtered,
            inflame_key,
            how='inner',
            predicate='intersects',
            target_coordinate_system=fullres_coord_key
        )
        if not joined_inflamed.empty:
            ids_to_discard = joined_inflamed.index.unique()
            mask_discard = tmp_sdata.tables[table_key].obs.index.isin(ids_to_discard)
            tmp_sdata.tables[table_key].obs.loc[mask_discard, 'to_discard'] = True
            print(f"Marked {mask_discard.sum()} nuclei as 'to_discard'")
        else:
            print("No nuclei found inside 'inflamed_poly'")
    else:
        print(f"Skipping 'to_discard':  '{inflame_key}' not found in shapes")
    # --- 4) Summary ---
    print("\nSummary:")
    print(f"  in_treatment: {tmp_sdata.tables[table_key].obs['in_treatment']. sum()} nuclei")
    print(f"  to_discard:    {tmp_sdata.tables[table_key].obs['to_discard'].sum()} nuclei")
    final_path = f"{save_path}{blocco_key}_{tissue_type}"
    tmp_sdata.write(final_path)

print("\nAll subsets processed and saved successfully.")
