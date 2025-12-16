import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import sopa
from skimage import io
from anndata import AnnData
import os

from py_scripts.utils.utils_fun import read_from_json
import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt

from PIL import Image
from spatialdata.transformations import Identity, Scale, set_transformation
from shapely.geometry import Polygon

"""
Pipeline vecchia:
    - create_zarr del blocco
    - bbox query e crop
    - ritaglio sugli intissue ecc
    
Pipeline nuova:
    - create_zarr del blocco
    - divisione nuclei per tissue_type
    - query per tissue_type e crop
    
"""

"""## **Helper Functions**
We use two custom helper function in this guide, which are defined in the next set of code cells.

The function `create_zarr` takes the raw output files from 10x Genomics' Visium HD processing and structures them into a single Zarr file, making the data ready for spatial analysis using libraries like `spatialdata`. It takes input paths, loads and prepares data, defines coordinate systems and transformations, processes the cell segmentation, integrates it with the `AnnData` object, creates `SpatialData` elements, and writes all of this to a Zarr file.

It takes five inputs:
* **`image_path`**: The path to the image file.
* **`count_matrix_path`**: The path to the cell segmentation filtered feature-barcode count matrix file.
* **`scale_factors_path`**: The path to the scale factors JSON file.
* **`geojson_path`**: The path to the cell segmentation GeoJSON file.
* **`intissue_geojson_path`**: The path to the intissue GeoJSON file.
* **`sample_name`**: A name for the Zarr output file.

In this guide, we use the `hires_image.png` file to reduce the size of the Zarr file and the RAM requirement, which speeds up the code. The full-resolution microscope image can also be used, with adjustments to the code to properly scale the cell segmentation GeoJSON file. The images required in the final Zarr file will depend on the downstream analyses.
"""


def create_zarr(count_matrix_path,
                image_path,
                scale_factors_path,
                geojson_path,
                intissue_geojson_path,
                sample_name,
                save_path
):
    print(sample_name)
    # Load and Prepare Raw Data
    # Define file paths
    COUNT_MATRIX_PATH = count_matrix_path
    IMAGE_PATH = image_path
    SCALE_FACTORS_PATH = scale_factors_path
    GEOJSON_PATH = geojson_path
    INTISSUE_GEOJSON_PATH = intissue_geojson_path
    # Load AnnData
    adata = sc.read_10x_h5(COUNT_MATRIX_PATH)
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    adata.obs.index = sample_name +"_" + adata.obs.index.astype(str)
    # Load and preprocess image data
    image_data = np.array(Image.open(IMAGE_PATH))
    if image_data.ndim == 2:
        image_data = image_data[np.newaxis, :, :] # Add channel dimension for grayscale
    elif image_data.ndim == 3:
        image_data = np.transpose(image_data, (2, 0, 1)) # (H, W, C) -> (C, H, W) for spatialdata
    # Load scale factors
    with open(SCALE_FACTORS_PATH, 'r') as f:
        scale_data = json.load(f)
    # Load GeoJSON data
    with open(GEOJSON_PATH, 'r') as f:
        geojson_data = json.load(f)
    # Load GeoJSON data
    intissue_geojson_data = gpd.read_file(INTISSUE_GEOJSON_PATH)
    # Define coordinate systems:
    # `downscale_to_hires`: The coordinate system where shapes are located, scaled relative to the hires resolution.
    hires_scale = scale_data['tissue_hires_scalef']
    # Transformation for shapes (from pixel to downscale_to_hires)
    shapes_transformations = {
       "downscale_to_hires": Scale(np.array([hires_scale, hires_scale]), axes=("x", "y")) # if the high-resolution microscope image is being used and Identity() transform would be performed.
    }
    # Transformation for the 'hires_tissue_image' (it's already in the 'downscale_to_hires' space visually)
    image_transformations = {
        "downscale_to_hires": Identity()
    }
    # Transformation for the 'intissue_geojson_data'
    geojson_transformations = {
        "downscale_to_hires": Identity()
    }
    # Process the polygons of the three tissues
    intissue_geojson_data = intissue_geojson_data.set_crs(None, allow_override=True)
    # Process Cell Segmentation (GeoJSON) and Integrate with AnnData
    # Create a mapping from adata.obs.index to geojson features
    geojson_features_map = {
        f"{sample_name}_cellid_{feature['properties']['cell_id']:09d}-1": feature
        for feature in geojson_data['features']
    }
    # Prepare data for GeoDataFrame and update adata.obs
    geometries = []
    cell_ids_ordered = []
    for obs_index_str in adata.obs.index:
        feature = geojson_features_map.get(obs_index_str)
        if feature:
            # Create shapely Polygon from coordinates
            polygon_coords = np.array(feature['geometry']['coordinates'][0])
            geometries.append(Polygon(polygon_coords))
            cell_ids_ordered.append(obs_index_str)
        else:
            geometries.append(None) # Or a suitable placeholder
            cell_ids_ordered.append(obs_index_str)
    # Remove None entries if any (or handle them upstream)
    valid_indices = [i for i, geom in enumerate(geometries) if geom is not None]
    geometries = [geometries[i] for i in valid_indices]
    cell_ids_ordered = [cell_ids_ordered[i] for i in valid_indices]
    # Create GeoDataFrame for shapes
    shapes_gdf = gpd.GeoDataFrame({
        'cell_id': cell_ids_ordered,
        'geometry': geometries
    }, index=cell_ids_ordered)
    # Update adata.obs with cluster information and spatial identifiers
    adata.obs['cell_id'] = adata.obs.index
    adata.obs['region'] = sample_name + '_nuclei_boundaries'
    adata.obs['region'] = adata.obs['region'].astype('category')
    adata = adata[shapes_gdf.index].copy() # Filter adata to match shapes_gdf
    # Define names for SpatialData elements
    IMAGE_KEY =  sample_name + '_hires_tissue_image'
    TABLE_KEY =  'segmentation_counts'
    SHAPES_KEY = sample_name + '_nuclei_boundaries'
    SHAPES_IT_KEY = sample_name + '_intissue_boundaries'
    # Create SpatialData elements directly
    sdata = spd.SpatialData(
        images={
            IMAGE_KEY: Image2DModel.parse(image_data, transformations=image_transformations)
        },
        tables={
            TABLE_KEY: TableModel.parse(
                adata,
                region=SHAPES_KEY, # Link table to shapes element
                region_key='region', # Column in adata.obs indicating region name
                instance_key='cell_id' # Column in adata.obs with instance IDs (cell_id)
            )
        },
        shapes={
            SHAPES_KEY: ShapesModel.parse(shapes_gdf, transformations=shapes_transformations),
            SHAPES_IT_KEY: ShapesModel.parse(intissue_geojson_data, transformations=geojson_transformations)
        }
    )
    
    # Add cell_type info to the nuclei shapes
    sdata[SHAPES_KEY] = sopa.spatial.sjoin(sdata,
      SHAPES_KEY,
      SHAPES_IT_KEY,
      how = "left",
      predicate = "intersects",
      target_coordinate_system = "downscale_to_hires"
    )
    tissue_mapping = sdata.shapes[SHAPES_KEY]['name']
    sdata.tables[TABLE_KEY].obs['tissue_type'] = sdata.tables[TABLE_KEY].obs.index.map(tissue_mapping)
   
    # filter out the NaN values AKA non inside tissue
    sdata.tables[TABLE_KEY].obs['in_tissue'] = sdata.tables[TABLE_KEY].obs['tissue_type'].notna()
    
    # 3. Filter the data to keep only the cells where 'in_tissue' is True.
    # First, get the list of cell IDs (the index) to keep.
    cells_to_keep = sdata.tables[TABLE_KEY].obs[sdata.tables[TABLE_KEY].obs['in_tissue']].index
    
    # Now, filter both the table and the shapes in the SpatialData object.
    # This ensures the entire object remains consistent.
    sdata.tables[TABLE_KEY] = sdata.tables[TABLE_KEY][cells_to_keep, :].copy()
    sdata.shapes[SHAPES_KEY] = sdata.shapes[SHAPES_KEY].loc[cells_to_keep].copy()
    
    # 4. As a best practice, convert the 'tissue_type' column to a category.
    # After filtering, there are no NaN values left, so this is safe and efficient.
    sdata.tables[TABLE_KEY].obs['tissue_type'] = sdata.tables[TABLE_KEY].obs['tissue_type'].astype('category')
    
    sdata.write(f"{save_path}{sample_name}", overwrite=True)
    del sdata
    gc.collect()

"""
The second helper function `crop0` ensures that the images generated from the analysis are cropped 
to the region of interest, aligning with the Visium HD Capture Area of each sample. 
It takes as input a `SpatialData` object (`x`), a target coordinate system (`crs`), 
and a bounding box dictionary (`bbox`). A bounding box dictionary is a way to represent 
a rectangular region in a 2D space using a Python dictionary. 
It typically contains the minimum and maximum coordinates for both the x and y axes that define 
the boundaries of the box. The function assumes that the `bbox` dictionary was created 
using `spd.get_extent` with the same coordinate system. 
Internally, the function calls `spd.bounding_box_query` and uses the minimum and maximum 
coordinates from the dictionary to subset the data from the `SpatialData` object 
that falls within this defined rectangle. This ensures that subsequent visualizations 
or analyses are focused only on the relevant part of the data. This is required because 
the microscope image is often larger than the Visium HD Gene Expression capture area.
"""

def crop0(x,crs,bbox):
    return spd.bounding_box_query(
        x,
        min_coordinate=[bbox['x'][0], bbox['y'][0]],
        max_coordinate=[bbox['x'][1], bbox['y'][1]],
        axes=("x", "y"),
        target_coordinate_system=crs,
    )


"""
Function to divide the 3 tissue from the blocco 
"""

def query_and_crop(sdata, save_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/",
                          poly_path = "/mnt/europa/valerio/data/json/geojson_dir/in_tissue/fullres/",
                          fullres_path = "/mnt/europa/valerio/HE_images/color_corrected/",
                          json_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell/"
):
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
    
    # Add the 'all_poly
    poly = gpd.read_file(f"{poly_path}{blocco_key}_allpolys_fullres.geojson")
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
    
    # Get all unique tissue types to iterate over
    all_tissue_types = sdata[table_key].obs['tissue_type'].unique().tolist()
    print(f"Found tissue types in '{blocco_key}': {all_tissue_types}")

    # --- 2. Loop through each tissue type and apply the script's logic ---
    for tissue_type in all_tissue_types:
        print(f"\n--- Processing tissue: {tissue_type} ---")

        # Define keys for the new, filtered elements for this iteration
        nuclei_key_filtered = f"{blocco_key}_{tissue_type}_nuclei_boundaries"
        intissue_key_filtered = f"{blocco_key}_{tissue_type}_intissue_boundaries"
        fullimage_key_filtered = f"{blocco_key}_{tissue_type}_full_image"
        poly_key = "various_poly"
        
        # a. Filter shapes and add them as NEW elements to the main sdata object
        sdata[intissue_key_filtered] = sdata[all_poly_key_orig][sdata[all_poly_key_orig]['name'] == tissue_type].copy()
        sdata[nuclei_key_filtered] = sdata[nuclei_key_orig][sdata[nuclei_key_orig]['name'] == tissue_type].copy()
        sdata[poly_key] = sdata[all_poly_key_orig][sdata[all_poly_key_orig]['name'].str.contains('fibre_trattate|infiammazione')]
        # b. Get the extent of the newly created filtered tissue boundary
        extent = spd.get_extent(
            sdata, 
            coordinate_system=blocco_key, 
            exact=True, 
            elements=[intissue_key_filtered]
        )

        # c. Crop the entire sdata object (including the newly added elements)
        cropped_sdata = crop0(sdata, blocco_key, extent)
        
        # d. Assemble a temporary sdata object from the cropped elements
        # Note: The table here is cropped but not yet filtered by tissue_type, matching your script
        tmp_sdata = spd.SpatialData(
            images=cropped_sdata.images,
            shapes={
                nuclei_key_filtered: cropped_sdata[nuclei_key_filtered],
                intissue_key_filtered: cropped_sdata[intissue_key_filtered],
                poly_key: cropped_sdata[poly_key]
            },
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
        centers = spd.get_centroids(tmp_sdata[nuclei_key_filtered], coordinate_system=coord_key).compute()
        tmp_sdata[table_key].obs['x'] = centers['x']
        tmp_sdata[table_key].obs['y'] = centers['y']
        tmp_sdata[table_key].obs['sample_id'] = f"{blocco_key}_{tissue_type}"

        # g. Save the final processed object
        
        print(f"Saving sample {blocco_key}_{tissue_type} to: {save_path}")
        final_path = f"{save_path}{blocco_key}_{tissue_type}"
        tmp_sdata.write(final_path)

    print("\nAll subsets processed and saved successfully.")


def query_and_crop_better(sdata, 
                          save_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/try/",
                          poly_path = "/mnt/europa/valerio/data/json/geojson_dir/in_tissue/fullres/",
                          fullres_path = "/mnt/europa/valerio/HE_images/color_corrected/",
                          fluo_path = "/mnt/europa/valerio/Fluo_images/warped_tif/",
                          json_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell/"
):
    """
    Iteratively creates and saves cropped subsets of a SpatialData object,
    following the user's specified procedural script logic for each tissue type.

    Args:
        sdata (spd.SpatialData): The input SpatialData object.
        save_path (str): The base directory where the output .zarr files will be saved.
        poly_path (str): The base directory where the geojson files are stored
        json_path (str): The base directory where the scalefactors are stored
        fullres_path (str): The directory of the H&E images in full resolution is stored
        fluo_path (str): The directory of the aligned fluorescence image is stored.
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
        cropped_sdata = crop0(sdata, fullres_coord_key, extent)
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
                ids_in_treatment = joined_treatment. index.unique()
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

#
def save_table_as_h5ad(spatial_data_path, table_key, output_dir):
    """
    Reads a SpatialData object from a folder, extracts a table, and saves it as an .h5ad file
    with a filename derived from the `stem` of SpatialData's path.

    Parameters:
        spatial_data_path (str): Path to the folder containing the SpatialData Zarr store.
        table_key (str): The key name of the table to extract from SpatialData.
        output_dir (str): Directory to save the .h5ad file.

    Returns:
        None
    """
    # Step 1: Read the SpatialData object
    print(f"Loading SpatialData object from: {spatial_data_path}")
    sdata = spd.read_zarr(spatial_data_path)

    # Step 2: Check if the table_key exists in tables
    if table_key not in sdata.tables:
        raise KeyError(f"Table '{table_key}' not found in SpatialData. Available tables are: {list(sdata.tables.keys())}")

    # Step 3: Extract the table as AnnData
    print(f"Extracting table '{table_key}'...")
    adata = sdata.tables[table_key]

    # Step 4: Generate the output path
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Derive the file name using the stem of sdata.path and append `.h5ad`
    file_name = f"{sdata.path.stem}.h5ad"
    output_path = os.path.join(output_dir, file_name)

    # Step 5: Save as .h5ad
    print(f"Saving table '{table_key}' as .h5ad to: {output_path}")
    adata.write_h5ad(output_path)

    print("Save complete!")
