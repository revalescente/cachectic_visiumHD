# required packages
# pip install anndata==0.12.0 pydeseq2==0.5.2 squidpy==1.6.5 spatialdata[extra]==0.4.0 geosketch==1.3 harmonypy==0.0.10 igraph==0.11.8

import spatialdata as spd
import spatialdata_plot as splt
import spatialdata_io as so
import geosketch as sketch
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce

import json
import gc
import geopandas as gpd
from spatialdata.models import Image2DModel, TableModel, ShapesModel
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from PIL import Image
from spatialdata.transformations import Identity, Scale
from shapely.geometry import Polygon

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
                sample_name
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
    SHAPES_IT_KEY = sample_name + 'intissue_boundaries'
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
            SHAPES_KEY: ShapesModel.parse(shapes_gdf, transformations=shapes_transformations)
            SHAPES_IT_KEY: ShapesModel.parse(intissue_geojson_data, transformations=geojson_transformations)
        }
    )
    sdata.write(sample_name, overwrite=True)
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

"""# **Section 3: Conversion of Space Ranger Output to Zarr Format and SpatialData Object Creation**

We begin this section by converting the standard Visium HD output into a Zarr file, 
as the `spatialdata` library expects data in this format. If you are running this code locally, 
this code needs to be run only once, as the datasets can be loaded directly from the saved Zarr files afterward.

In the following code snippet, a dictionary is created where each sample name serves as a unique key. 
The value associated with each key is a list containing:
* the path and filename of the filtered feature-cell matrix in `h5` format.
* the location and name of the image to be stored in the `SpatialData` object.
* the scale factors JSON file, so that the cell segmentation results can beoverlaid onto the tissue image.
* the cell segmentation GeoJSON file, so the cell segmentation results can be visualized on the tissue image.
* the desired name for the Zarr file.

Each key-value pair in this dictionary is then processed by the `create_zarr` helper function. 
For this specific example, we only use the `tissue_hires_image.png`. 
Other images, such as a high-resolution microscope image or CytAssist image, can be added to the `SpatialData` object.

"""
# /mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs
# /mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE

# Create and save Zarr files for the cell segmentation outputs.
sample = {
  "blocco1":["/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/raw_feature_cell_matrix.h5",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/spatial/tissue_hires_image.png",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/spatial/scalefactors_json.json",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/nucleus_segmentations.geojson",
             "/mnt/europa/valerio/data/json/geojson_dir/tissue_hires_image_blocco1.geojson",
             "blocco1_nuclei"]
}

print("Saving zarr files")
for key, inputs in sample.items():
    create_zarr(count_matrix_path=inputs[0],
                image_path=inputs[1],
                scale_factors_path=inputs[2],
                geojson_path=inputs[3],
                intissue_geojson_path=inputs[4],
                sample_name=inputs[5])

del sample, inputs, key
gc.collect()

"""
## **The SpatialData Object and Its Components**
The `spatialdata` library is built to manage and analyze multiomic spatial datasets. It brings together multiple 
data types into a single, unified `SpatialData` object. These objects act as on-disk containers that utilize the 
Zarr file format to store various **Elements** or data types.
A `SpatialData` object created from Visium HD data typically has the following Elements:
* **Images**: CytAssist and microscopy images (e.g., H&E, fluorescence) providing spatial context. 
These can be accessed via `sdata.images`.
* **Shapes**: Geometric annotations such as polygons or circles representing regions of interest, cells, or spots. 
In Visium HD, these often represent binned regions or cell segmentations. These are accessible via `sdata.shapes`.
* **Tables**: An `AnnData` object associated with the spatial elements, typically containing gene expression data, 
cellular metadata, and computational results (e.g., clusters, UMAP embeddings). This Element is used for 
downstream analyses and is accessed via `sdata.tables`. Each `AnnData` table within `sdata.tables` has:
  * `.X`: The primary data matrix (e.g., raw counts, normalized counts).
  * `.obs`: Observation metadata (e.g., sample ID, cluster assignments).
  * `.var`: Variable metadata (e.g., gene names).
  * `.obsm`: Multi-dimensional annotations (e.g., PCA, UMAP embeddings).
  * `.layers`: Alternative representations of `.X`.

In the next code block, the `SpatialData` object is created from the Zarr files. 
Each Zarr file is read into a list of `SpatialData` objects using the `read_zarr` 
function before the objects are concatenated. In addition, a `sample` column is added to 
each `AnnData` table within the `SpatialData` object, and the `var_names_make_unique` function 
is used to ensure that gene names are unique.
"""

sdata = spd.read_zarr("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/blocco1_nuclei")

# divide the 3 tissue

s1 = sdata2.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 0],
    max_coordinate=[4185, 1800],
    target_coordinate_system="downscale_to_hires",
)
s1_extent = spd.get_extent(s1, coordinate_system='downscale_to_hires', exact=True, 
                  has_images=True, 
                  has_labels=True, 
                  has_points=True, 
                  has_shapes=True, 
                  elements=['blocco1_nuclei_intissue_boundaries']
)
s1_crop = crop0(s1, "downscale_to_hires", s1_extent)
#
s2 = sdata2.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 2000],
    max_coordinate=[4185, 4000],
    target_coordinate_system="downscale_to_hires",
)
s2_extent = spd.get_extent(s2, coordinate_system='downscale_to_hires', exact=True, 
                  has_images=True, 
                  has_labels=True, 
                  has_points=True, 
                  has_shapes=True, 
                  elements=['blocco1_nuclei_intissue_boundaries']
)

s2_crop = crop0(s2, "downscale_to_hires", s2_extent)
#
s3 = sdata2.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[0, 4000],
    max_coordinate=[4185, 6000],
    target_coordinate_system="downscale_to_hires",
)
s3_extent = spd.get_extent(s3, coordinate_system='downscale_to_hires', exact=True, 
                  has_images=True, 
                  has_labels=True, 
                  has_points=True, 
                  has_shapes=True, 
                  elements=['blocco1_nuclei_intissue_boundaries']
)

s3_crop = crop0(s3, "downscale_to_hires", s3_extent)


plt.figure(figsize=(50, 50))
ax = plt.gca()
s1_crop.pl.render_images("blocco1_nuclei_hires_tissue_image").pl.render_shapes("blocco1_nuclei_nuclei_boundaries", outline=True, outline_alpha=1, outline_width=3, fill_alpha=0
).pl.show(ax = ax, coordinate_systems = "downscale_to_hires", save = "output_python/sranger/blocco1_s1.png")


s1_crop.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_c26STAT3")
s2_crop.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_sham")
s3_crop.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_c26foxO")

# saved but i need to check the intissue 

tissue_type = "c26STAT3"
image_key = "blocco1_nuclei_hires_tissue_image"
intissue_key = "blocco1_nuclei_intissue_boundaries"
nuclei_key_old = "blocco1_nuclei_nuclei_boundaries"
nuclei_key = f"blocco1_{tissue_type}_nuclei_boundaries"
table_key = "segmentation_counts"
coord_key = "downscale_to_hires"

s2_crop[intissue_key] = s2_crop[intissue_key].iloc[[1]]

# filter intissue
nuclei_intissue = sopa.spatial.sjoin(s2_crop,
      nuclei_key_old,
      intissue_key,
      how = "inner",
      predicate = "intersects",
      target_coordinate_system = coord_key
)

# remove duplicated
nuclei_intissue = nuclei_intissue[nuclei_intissue['cell_id'].duplicated(keep=False) == False]
nuclei_intissue = nuclei_intissue[['cell_id', 'geometry', 'name']]

s2_crop[nuclei_key] = nuclei_intissue

# Annotate the table with the shapes
s2_crop[table_key].obs["region"] = pd.Categorical([nuclei_key]*len(s2_crop[table_key]))
s2_crop[table_key].uns["spatialdata_attrs"] = {
    "region": nuclei_key,  # name of the Shapes element we will use later
    "region_key": "region",      # column in adata.obs that will link a given obs to the elements it annotates
    "instance_key": "cell_id",  # column that matches a given obs in the table to a given circle
}
s2_crop.set_table_annotates_spatialelement(table_key, region=nuclei_key)    
    
# Filter the table
s2_crop[table_key] = spd.match_table_to_element(s2_crop, 
    element_name=nuclei_key, 
    table_name=table_key
)

# Map the 'exp_condition' in the 'filtered' table
location_to_name = s2_crop[nuclei_key].set_index('cell_id')['name']
# Map the names to adata.obs using location_id
s2_crop[table_key].obs['tissue_type'] = s2_crop[table_key].obs['cell_id'].map(location_to_name)

# add nuclei centroids to .obs 
centers = spd.get_centroids(blocco1_c26foxO['blocco1_c26foxO_nuclei_boundaries'], coordinate_system='downscale_to_hires').compute()

blocco1_c26foxO[table_key].obs['x'] = centers['x']
blocco1_c26foxO[table_key].obs['y'] = centers['y']

# other useful things
blocco1_c26foxO[table_key].obs['intissue'] = True
blocco1_c26foxO[table_key].obs['sample_id'] = f"blocco1_{tissue_type}"

temp_tab = blocco1_c26foxO[table_key].copy()
blocco1_c26foxO.delete_element_from_disk(table_key)
del blocco1_c26foxO[table_key]
blocco1_c26foxO[table_key] = temp_tab
blocco1_c26foxO.write_element(table_key)
# save the three samples

blocco1_c26STAT3 = s1_crop.subset(element_names = [image_key, nuclei_key, intissue_key, table_key])
blocco1_sham = s2_crop.subset(element_names = [image_key, nuclei_key, intissue_key, table_key])
blocco1_c26foxO = s3_crop.subset(element_names = [image_key, nuclei_key, intissue_key, table_key])

blocco1_c26STAT3.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_c26STAT3")
blocco1_sham.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_sham")
blocco1_c26foxO.write("/mnt/europa/valerio/data/zarr_store/spaceranger_v4/samples/blocco1_c26foxO")
