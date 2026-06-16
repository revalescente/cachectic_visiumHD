import geopandas as gpd
import numpy as np
import rasterio.features
import shapely.geometry


def mask_to_polygons(mask):
    geometries = []
    labels = []

    for geometry, value in rasterio.features.shapes(
        mask.astype(np.int32),
        mask=mask > 0,
        connectivity=8,
    ):
        geometries.append(shapely.geometry.shape(geometry))
        labels.append(int(value))

    polygons = gpd.GeoDataFrame(
        {"geometry": geometries, "label_id": labels},
        index=labels,
    )

    if len(labels) != len(set(labels)):
        polygons = polygons.dissolve(by="label_id")
        polygons["label_id"] = polygons.index

    return polygons
