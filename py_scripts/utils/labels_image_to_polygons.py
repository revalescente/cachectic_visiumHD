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

# convert to tiff from ome.tiff
# ./bftools/bfconvert -series 0 -bigtiff \
#   data/manual_annotations/train_images_100_regions_background.ome.tiff \
#   data/manual_annotations/tiff/train_images_100_regions_background.tiff

# read the manual annotation - nuclei
label_path = "/mnt/europa/valerio/data/manual_annotations/tiff/train_images_100_regions_selezioni.tiff"
selezioni = skimage.io.imread(label_path).astype(int)

print(f"Number of unique objects: {len(np.unique(selezioni)) - 1}") # subtract 1 for background (0)

# read the composition annotated image
img_path = "/mnt/europa/valerio/data/manual_annotations/train_images_100_regions.tif"
img = skimage.io.imread(img_path)

# create the sd object with the image
img_parsed = Image2DModel.parse(data=img, 
            scale_factors=(2, 2, 2), 
            transformations={"standard": Identity()},
            dims=("y", "x", "c")
)  

selezioni_parsed = Labels2DModel.parse(data=selezioni, 
            transformations={"standard": Identity()},
            dims=("y", "x")
)  

train_sd = sd.SpatialData(images={"fluo_image": img_parsed}, labels={"nuclei_annots": selezioni_parsed})

# plottino

plt.figure(figsize=(20, 20))
ax = plt.gca()
train_sd.pl.render_labels('nuclei_annots'
).pl.show(ax = ax, coordinate_systems="standard", save = 'output_python/manual_annots_labels.png')

# transformation to polygons

train_sd["nuclei_annots_poly"] = to_polygons(train_sd["nuclei_annots"])

# extract from sd in the correct coord system
selezioni_poly = train_sd.transform_element_to_coordinate_system(
    "nuclei_annots_poly", "standard"
)

# 3. Export to GeoJSON
output_path = "/mnt/europa/valerio/data/manual_annotations/geojson/nuclei_annotations_poly.geojson"
selezioni_poly.to_file(output_path, driver='GeoJSON', index=False)

# ------

# read the manual annotation - nuclei
label_path = "/mnt/europa/valerio/data/manual_annotations/tiff/train_images_100_regions_background.tiff"
background = skimage.io.imread(label_path).astype(int)

print(f"Number of unique objects: {len(np.unique(selezioni)) - 1}") # subtract 1 for background (0)



background_parsed = Labels2DModel.parse(data=background, 
            transformations={"standard": Identity()},
            dims=("y", "x")
)  

train_sd["background_annots"] = background_parsed

# transformation to polygons

train_sd["background_annots_poly"] = to_polygons(train_sd["background_annots"])

# extract from sd in the correct coord system
background_poly = train_sd.transform_element_to_coordinate_system(
    "background_annots_poly", "standard"
)

# 3. Export to GeoJSON
output_path = "/mnt/europa/valerio/data/manual_annotations/geojson/background_annotations_poly.geojson"
background_poly.to_file(output_path, driver='GeoJSON', index=False)
