import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import geopandas as gpd
import pandas as pd
import skimage
import numpy as np
from spatialdata.models import ShapesModel, Image2DModel, Labels2DModel
from spatialdata import to_polygons
from spatialdata.transformations import Identity, set_transformation, get_transformation
from py_scripts.utils.utils_fun import read_from_json
import py_scripts.segmentation.segm_functions as sf
from skimage.measure import regionprops_table

# convert to tiff from ome.tiff
# ./bftools/bfconvert -series 0 -bigtiff \
#   data/manual_annotations/train_images_100_regions_background.ome.tiff \
#   data/manual_annotations/tiff/train_images_100_regions_background.tiff

# read the manual annotation - nuclei
label_path = "/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks/blocco1_sham_finalprediction.tif"
nuclei_arivis = skimage.io.imread(label_path).astype(int)

print(f"Number of unique objects: {len(np.unique(selezioni)) - 1}") # subtract 1 for background (0)

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco1_sham.zarr")

# i need to transform the nuclei into the full_image coordinates system so i need to extract it first

transf = get_transformation(sdata["full_image"], "blocco1")

nuclei_arivis_parsed = Labels2DModel.parse(data=nuclei_arivis, 
            transformations={"blocco1": transf},
            dims=("y", "x")
)  
sdata["label_nuclei_arivis"] = nuclei_arivis_parsed

# plottino

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[8000, 13000],
    max_coordinate=[9000, 14000],
    target_coordinate_system="blocco1",
).pl.render_images('full_image'
).pl.render_labels('label_nuclei_arivis', color = "green"
).pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/arivis_nuclei.png')

# everything in the right place, it seems

# transformation to polygons

sdata["nuclei_arivis_poly"] = to_polygons(sdata["label_nuclei_arivis"])

# aggregate to obtain the matrix of gene vs arivis_nuclei
sopa.utils.set_sopa_attrs(sdata, 
              cell_segmentation_key='full_image', 
              tissue_segmentation_key='full_image', 
              transcripts_key=None, 
              boundaries_key='nuclei_arivis_poly', 
              bins_table_key='filtered'
 )

sopa.aggregate(sdata, key_added = 'arivis_nuclei', bins_key= "filtered",
  shapes_key = "nuclei_arivis_poly", expand_radius_ratio=0, min_transcripts=10,
  min_intensity_ratio=0.15, no_overlap = True
)


element_extent = sd.get_extent(sdata['nuclei_arivis_poly'], coordinate_system='blocco1', exact=True)
sdata['raster_arivis_nuclei'] = sd.rasterize(
    sdata['nuclei_arivis_poly'],
    axes=["x", "y"],
    min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
    max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
    target_coordinate_system='blocco1',
    target_unit_to_pixels=1,
)
# 4. --- Prepare Label Mask ---
# Squeeze the array to (H, W) and ensure it's an integer type
label_mask = sdata['raster_arivis_nuclei'].values.squeeze().astype(np.int32)

properties_to_extract = [
    'label', 'area', 'eccentricity', 'solidity', 'extent',
    'major_axis_length', 'minor_axis_length'
]

props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))

label_to_id = sdata['raster_arivis_nuclei'].attrs['label_index_to_category']
props_df['cell_id'] = props_df['label'].map(label_to_id)
props_df = props_df.set_index('cell_id')
props_df = props_df.drop(columns='label')

cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
sdata['arivis_nuclei'].obs = sdata['arivis_nuclei'].obs.join(props_df[cols], how='left')

sdata.write_element('label_nuclei_arivis')
sdata.write_element('raster_arivis_nuclei')
sdata.delete_element_from_disk('arivis_nuclei')
sdata.write_element('arivis_nuclei')

#export adata
adata = sdata[table_name].copy()
adata.obs['x_coord'] = adata.obsm['spatial'][:, 0]
adata.obs['y_coord'] = adata.obsm['spatial'][:, 1]

del adata.obsm
adata.write_h5ad('/mnt/europa/valerio/data/zarr_store/adatas/arivis_segmentation_tables/blocco1_sham.h5ad')
