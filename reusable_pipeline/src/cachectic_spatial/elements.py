from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel

from cachectic_spatial.polygons import mask_to_polygons


def add_annotations(sdata, annotations, key, coordinate_system, transformation):
    sdata.shapes[key] = ShapesModel.parse(
        annotations,
        transformations={coordinate_system: transformation},
    )


def add_fluorescence_image(sdata, image, key, coordinate_system, transformation):
    sdata.images[key] = Image2DModel.parse(
        data=image,
        scale_factors=(2, 2, 2),
        transformations={coordinate_system: transformation},
        dims=("y", "x", "c"),
    )


def add_segmented_shapes(
    sdata,
    mask,
    labels_key,
    shapes_key,
    coordinate_system,
    transformation,
):
    sdata.labels[labels_key] = Labels2DModel.parse(
        data=mask,
        transformations={coordinate_system: transformation},
        dims=("y", "x"),
    )

    polygons = mask_to_polygons(mask)
    sdata.shapes[shapes_key] = ShapesModel.parse(
        polygons,
        transformations={coordinate_system: transformation},
    )
