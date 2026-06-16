import json

import geopandas as gpd
from skimage import io
from spatialdata_io import visium_hd

from cachectic_spatial.config import Paths, SegmentedObject


def read_samples(path):
    with open(path) as file:
        return json.load(file)


def read_visium(block, paths: Paths):
    return visium_hd(
        path=paths.visium.format(block=block),
        fullres_image_file=paths.fullres_image.format(block=block),
        dataset_id=block,
        filtered_counts_file=False,
        bin_size=["002", "008", "016"],
        bins_as_squares=True,
        annotate_table_by_labels=False,
        load_all_images=False,
        var_names_make_unique=True,
        image_models_kwargs={"dims": ["c", "y", "x"]},
        load_segmentations_only=False,
        load_nucleus_segmentations=False,
    )


def read_annotations(block, sample, paths: Paths):
    annotations = gpd.read_file(paths.annotations.format(block=block, sample=sample))
    return annotations.set_crs(None, allow_override=True)


def read_fluorescence(block, sample, paths: Paths):
    image = io.imread(paths.fluorescence_image.format(block=block, sample=sample))
    return image[:, :, 0:1]


def read_segmentation(block, sample, object_spec: SegmentedObject):
    path = object_spec.segmentation_path.format(block=block, sample=sample)
    return io.imread(path).astype(int)
