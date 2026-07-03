import geopandas as gpd
import numpy as np
import rasterio.features
import shapely.geometry
from skimage.color import rgb2gray
from skimage.filters import threshold_otsu
from skimage.morphology import closing, disk, remove_small_holes, remove_small_objects


def background_mask(image):
    gray = rgb2gray(image)
    threshold = threshold_otsu(gray)

    tissue = gray < threshold
    tissue &= image.max(axis=2) > 20
    tissue = closing(tissue, disk(3))
    tissue = remove_small_objects(tissue, max_size=100)
    tissue = closing(tissue, disk(3))
    tissue = remove_small_holes(tissue, max_size=10000)

    return ~tissue


def mask_to_polygons(mask, label="background"):
    geometries = []
    values = []
    mask = np.asarray(mask).astype(np.uint8)
    for geometry, value in rasterio.features.shapes(
        mask,
        mask=mask > 0,
        connectivity=8,
    ):
        geometries.append(shapely.geometry.shape(geometry))
        values.append(label if value else None)

    gdf = gpd.GeoDataFrame({"label": values, "geometry": geometries})
    gdf = gdf[gdf["label"].notna()].reset_index(drop=True)
    return gdf.set_crs(None, allow_override=True)


def save_background_geojson(mask, output):
    gdf = mask_to_polygons(mask)
    gdf.to_file(output, driver="GeoJSON", index=False)
    return gdf
