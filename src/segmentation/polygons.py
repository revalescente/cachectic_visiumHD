from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from rasterio import features


def to_polygons(mask):
    mask = mask.astype(np.int32)

    geometries = []
    label_ids = []

    for geometry, label_id in features.shapes(
        mask,
        mask=mask > 0,
        connectivity=8,
    ):
        geometries.append(shapely.geometry.shape(geometry))
        label_ids.append(int(label_id))

    gdf = gpd.GeoDataFrame(
        {"label_id": label_ids, "geometry": geometries},
        index=pd.Index(label_ids),
    )

    gdf = gdf.dissolve(by="label_id", as_index=True)
    gdf["label_id"] = gdf.index

    return gdf


def save_gdf(gdf, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(output_path, driver="GeoJSON")
