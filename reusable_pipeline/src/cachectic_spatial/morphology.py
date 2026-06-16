import numpy as np
import pandas as pd
import spatialdata as sd
from skimage.measure import regionprops_table


def extract_morphology(
    sdata,
    shapes_key,
    raster_key,
    coordinate_system,
    properties,
):
    extent = sd.get_extent(
        sdata[shapes_key],
        coordinate_system=coordinate_system,
        exact=True,
    )

    sdata[raster_key] = sd.rasterize(
        sdata[shapes_key],
        axes=["x", "y"],
        min_coordinate=[extent["x"][0], extent["y"][0]],
        max_coordinate=[extent["x"][1], extent["y"][1]],
        target_coordinate_system=coordinate_system,
        target_unit_to_pixels=1,
    )

    label_mask = sdata[raster_key].values.squeeze().astype(np.int32)
    data = pd.DataFrame(
        regionprops_table(label_mask, properties=("label", *properties))
    )

    label_to_id = sdata[raster_key].attrs["label_index_to_category"]
    data["cell_id"] = data["label"].map(label_to_id)

    return data.set_index("cell_id").drop(columns="label")


def add_morphology_to_table(sdata, table_key, morphology):
    sdata.tables[table_key].obs = sdata.tables[table_key].obs.join(
        morphology,
        how="left",
    )
