import spatialdata as sd
from spatialdata_io import visium_hd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import re
import os
import numpy as np
import geopandas as gpd
import pandas as pd
from skimage import io
from shapely.affinity import scale
from spatialdata.models import ShapesModel, Image2DModel, Labels2DModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from skimage.measure import regionprops_table
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.pp_sdata.pp_functions as pp
import py_scripts.segmentation.segm_functions as sf

# disable autosave of sopa
sopa.settings.auto_save_on_disk = False

# keys name
BLOCCO_KEY = 'blocco3'
# SAMPLE_KEY = 'c26murf1'
#HIRES_KEY = f'{BLOCCO_KEY}_hires_image'
#FULL_IMAGE_KEY = 'full_image'
#FLUO_KEY = 'fluo_image'

# paths of interest
intissue_gfp_dir = "/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys"
data_blocco2 = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell"
data_b_all = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_cellexpans"
arivis_dir = "/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks"
fullres_dir = "/mnt/europa/valerio/HE_images/color_corrected/blocchi"
fluo_path = "/mnt/europa/valerio/Fluo_images/warped_tif/samples"
save_dir = "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins"

# full res image path
fullres_path = f"{fullres_dir}/pp_{BLOCCO_KEY}_20x.tif"

# 1. load the data
sdata = visium_hd(
        path = f"{data_b_all}/{BLOCCO_KEY}/outs",
        fullres_image_file = fullres_path,
        dataset_id = BLOCCO_KEY,
        filtered_counts_file = False,
        bin_size = ['002','008','016'],
        bins_as_squares = True,
        annotate_table_by_labels = False,
        load_all_images = False,
        var_names_make_unique = True,
        image_models_kwargs = {'dims' : ['c', 'y', 'x']},
        load_segmentations_only = False,
        load_nucleus_segmentations = False
)

# 1b. remove the cell segmentations of spaceranger
if f'{BLOCCO_KEY}_cell_segmentations' in sdata:
    del sdata[f'{BLOCCO_KEY}_cell_segmentations']
if 'cell_segmentations' in sdata:
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
        # save path and check if already exist, in this case skip the sample
        save_path = f"{save_dir}/{BLOCCO_KEY}_{sample_name}"
        if os.path.isdir(save_path):
            print(f"  [SKIP] Output directory already exists for sample '{sample_name}': {save_path}")
            continue
        # 2. Extract the right transformation from the full_image
        transf = get_transformation(sdata_bbox[f"{BLOCCO_KEY}_full_image"], to_coordinate_system = BLOCCO_KEY)
        # 2a. Add intissue_GPF_polys to the sdata_bbox
        all_poly = gpd.read_file(intissue_gfp_dir + f"/{BLOCCO_KEY}_{sample_name}.geojson")
        all_poly = all_poly.set_crs(None, allow_override=True)
        # Parsed polys inside the sdata_bbox shapes
        polys_name = f"intissue_GFP_poly_{sample_name}"
        polys_parse = ShapesModel.parse(all_poly, transformations={BLOCCO_KEY: transf})
        sdata_bbox.shapes[polys_name] = polys_parse
        # 2b. Add the fluo image to the sdata_bbox
        fluo_image = io.imread(f"{fluo_path}/{BLOCCO_KEY}_{sample_name}_fluo_image.tif")
        fluo_image = fluo_image[:, :, 0:1] # shape will be (y, x, 1)
        # this creates 4 downscaled images as well
        fluo_parsed = Image2DModel.parse(data=fluo_image, 
                    scale_factors=(2, 2, 2), 
                    transformations={BLOCCO_KEY: transf},
                    dims=("y", "x", "c")
        )  
        sdata_bbox[f'{BLOCCO_KEY}_fluo_image'] = fluo_parsed
        # 3. Spatial join of the intissue and GFP polygons for the 3 tables
        table_shape_mapping = [
            ('square_016um', f'{BLOCCO_KEY}_square_016um'),
            ('square_008um', f'{BLOCCO_KEY}_square_008um'),
            ('square_002um', f'{BLOCCO_KEY}_square_002um')
        ]
        for table_name, shape_key in table_shape_mapping:
            print(f"  Processing table: {table_name}")
            left_element = sdata_bbox.transform_element_to_coordinate_system(shape_key, BLOCCO_KEY)
            right_element = sdata_bbox.transform_element_to_coordinate_system(polys_name, BLOCCO_KEY)
            joined_bins = gpd.sjoin(left_element, right_element, how='inner', predicate='intersects')
            
            if table_name not in sdata_bbox.tables:
                print(f"    [ERROR] Table '{table_name}' not found.")
                continue
            table = sdata_bbox.tables[table_name]
            instance_key = table.uns['spatialdata_attrs']['instance_key']
            for col in ['in_treatment', 'to_discard', 'in_tissue']:
                table.obs[col] = False
            nuclei_per_poly = joined_bins.groupby('name').apply(lambda x: x.index.tolist())
            for source_key in nuclei_per_poly.index:
                col_name = None
                if "fibre_trattate" in source_key: col_name = "in_treatment"
                elif "infiammazione" in source_key: col_name = "to_discard"
                elif sample_name in source_key: col_name = "in_tissue"
                if col_name:
                    nuclei_list = nuclei_per_poly[source_key]
                    mask = table.obs[instance_key].isin(nuclei_list)
                    if mask.any():
                        table.obs.loc[mask, col_name] = True
                        print(f"    Mapped '{source_key}' -> '{col_name}' ({mask.sum()})")
            # Filter intissue bins
            filtered_table = table[table.obs['in_tissue'] == True].copy()
            sdata_bbox.tables[table_name] = filtered_table
            # matching shapes with filtered tables to align geometries
            matched_elements = sd.match_element_to_table(sdata_bbox, shape_key, table_name)
            sdata_bbox.shapes[shape_key] = matched_elements[0][shape_key]
            # 4. Calculate GFP values 
            channel_aggregation = sopa.aggregation.aggregate_channels(
                sdata_bbox, 
                image_key=f'{BLOCCO_KEY}_fluo_image', 
                shapes_key=shape_key, 
                expand_radius_ratio=0, 
                mode='max', 
                no_overlap=False
            ) 
            max_values_vector = channel_aggregation.max(axis=1)
            # Add to the table, handling if we only want it for `in_treatment`
            # (By default this adds it to all filtered bins; you can mask it later if needed)
            sdata_bbox.tables[table_name].obs['GFP_value'] = max_values_vector
            print(f"  Finished processing {table_name}")
        # 5. ADD ARIVIS NUCLEI
        print(f"  Adding Arivis nuclei for {sample_name}...")
        arivis_path = f"{arivis_dir}/{BLOCCO_KEY}_{sample_name}_finalprediction.tiff"
        nuclei_arivis = io.imread(arivis_path).astype(int)
        num_objects = len(np.unique(nuclei_arivis)) - 1 
        if num_objects < 0: 
            num_objects = 0
        print(f"    Number of arivis objects: {num_objects}")
        nuclei_arivis_parsed = Labels2DModel.parse(
            data=nuclei_arivis, 
            transformations={BLOCCO_KEY: transf},
            dims=("y", "x")
        )  
        sdata_bbox["label_nuclei_arivis"] = nuclei_arivis_parsed 
        # transformation to polygons
        nuclei_shapes = sf.precise_to_polygons(sdata_bbox["label_nuclei_arivis"])
        nuclei_shapes_parsed = ShapesModel.parse(nuclei_shapes, transformations={BLOCCO_KEY : transf})
        sdata_bbox["nuclei_arivis_poly"] = nuclei_shapes_parsed
        # aggregate to obtain the matrix of gene vs arivis_nuclei
        sopa.utils.set_sopa_attrs(
            sdata_bbox, 
            cell_segmentation_key=f'{BLOCCO_KEY}_full_image', 
            tissue_segmentation_key=f'{BLOCCO_KEY}_full_image', 
            transcripts_key=None, 
            boundaries_key='nuclei_arivis_poly', 
            bins_table_key='square_002um'
        )
        sopa.aggregate(
            sdata_bbox, 
            key_added='arivis_nuclei_table', 
            bins_key="square_002um",
            shapes_key="nuclei_arivis_poly", 
            expand_radius_ratio=0, 
            min_transcripts=10,
            min_intensity_ratio=0.15, 
            no_overlap=True
        )
        # GFP values for the arivis nuclei
        channel_aggregation = sopa.aggregation.aggregate_channels(
            sdata_bbox, 
            image_key=f'{BLOCCO_KEY}_fluo_image', 
            shapes_key='nuclei_arivis_poly', 
            expand_radius_ratio=0, 
            mode='max', 
            no_overlap=False
        )
        max_values_vector = channel_aggregation.max(axis=1)
        sdata_bbox.tables['arivis_nuclei_table'].obs['GFP_value'] = max_values_vector
        # 6. FEATURE EXTRACTION
        print("  Extracting region properties...")
        element_extent = sd.get_extent(sdata_bbox['nuclei_arivis_poly'], coordinate_system=BLOCCO_KEY, exact=True)
        sdata_bbox['raster_arivis_nuclei'] = sd.rasterize(
            sdata_bbox['nuclei_arivis_poly'],
            axes=["x", "y"],
            min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
            max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
            target_coordinate_system=BLOCCO_KEY,
            target_unit_to_pixels=1,
        )
        label_mask = sdata_bbox['raster_arivis_nuclei'].values.squeeze().astype(np.int32)
        properties_to_extract = [
            'label', 'area', 'eccentricity', 'solidity', 'extent',
            'major_axis_length', 'minor_axis_length'
        ]
        props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))
        label_to_id = sdata_bbox['raster_arivis_nuclei'].attrs['label_index_to_category']
        props_df['cell_id'] = props_df['label'].map(label_to_id)
        props_df = props_df.set_index('cell_id')
        props_df = props_df.drop(columns='label')
        cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
        sdata_bbox.tables['arivis_nuclei_table'].obs = sdata_bbox.tables['arivis_nuclei_table'].obs.join(props_df[cols], how='left')
        # 7. SAVE
        sdata_bbox.write(save_path)
        print(f"  Successfully processed and saved {sample_name} to {save_path}\n")
