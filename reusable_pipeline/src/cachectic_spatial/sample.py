from spatialdata.transformations import get_transformation

from cachectic_spatial.annotations import annotate_table
from cachectic_spatial.config import BIN_SIZES, PATHS, SEGMENTED_OBJECTS
from cachectic_spatial.elements import add_annotations, add_fluorescence_image
from cachectic_spatial.fluorescence import add_fluorescence_to_table
from cachectic_spatial.input import read_annotations, read_fluorescence
from cachectic_spatial.objects import add_segmented_objects
from cachectic_spatial.output import write_sample


def process_sample(sdata, block, sample, sample_info):
    sample_data = sdata.query.bounding_box(
        axes=["x", "y"],
        min_coordinate=sample_info["min_coordinate"],
        max_coordinate=sample_info["max_coordinate"],
        target_coordinate_system=block,
    )

    image_key = f"{block}_full_image"
    transformation = get_transformation(sample_data[image_key], to_coordinate_system=block)

    annotations_key = "sample_annotations"
    annotations = read_annotations(block, sample, PATHS)
    add_annotations(sample_data, annotations, annotations_key, block, transformation)

    fluorescence_key = f"{block}_fluorescence"
    fluorescence = read_fluorescence(block, sample, PATHS)
    add_fluorescence_image(
        sample_data,
        fluorescence,
        fluorescence_key,
        block,
        transformation,
    )

    for bin_size in BIN_SIZES:
        table_key = f"square_{bin_size}um"
        shapes_key = f"{block}_square_{bin_size}um"
        annotate_table(
            sample_data,
            table_key,
            shapes_key,
            annotations_key,
            block,
            sample,
            filter_in_tissue=True,
        )
        add_fluorescence_to_table(
            sample_data,
            fluorescence_key,
            shapes_key,
            table_key,
        )

    for object_spec in SEGMENTED_OBJECTS:
        add_segmented_objects(
            sample_data,
            block,
            sample,
            object_spec,
            transformation,
            fluorescence_key,
            annotations_key,
        )

    write_sample(sample_data, block, sample, PATHS.output)
