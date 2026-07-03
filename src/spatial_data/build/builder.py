import os
import shutil
from pathlib import Path

import geopandas as gpd
from dask_image.imread import imread as dask_imread
from spatialdata.models import Image2DModel
from spatialdata.models import ShapesModel
from spatialdata.transformations import (
    Sequence,
    Translation,
    get_transformation,
    set_transformation,
)
from spatialdata_io import visium_hd

from src.geojson_parser import load_samples, rotated_rectangle, sample_bounds


def write_sdata(sdata, path, overwrite=True):
    path = Path(path)
    if path.exists() and overwrite:
        shutil.rmtree(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sdata.write(path)
    return path


def find_fullres_image_key(sdata, block):
    preferred = f"{block}_full_image"
    if preferred in sdata.images:
        return preferred
    fullres_like = [key for key in sdata.images if "full" in key.lower()]
    if fullres_like:
        return fullres_like[0]
    return next(iter(sdata.images))


def add_geojson_shapes(sdata, path, key, coordinate_system, transformation):
    os.environ.setdefault("OGR_GEOJSON_MAX_OBJ_SIZE", "0")
    gdf = gpd.read_file(path).set_crs(None, allow_override=True)
    gdf.attrs = {}
    sdata.shapes[key] = ShapesModel.parse(
        gdf,
        transformations={coordinate_system: transformation},
    )
    return gdf


def read_image(path):
    image = dask_imread(str(path))
    if image.ndim >= 3 and image.shape[0] == 1:
        image = image[0]
    return image


def image_dims(image):
    if image.ndim == 2:
        return ("y", "x")
    if image.ndim == 3 and image.shape[-1] in (1, 3, 4):
        return ("y", "x", "c")
    if image.ndim == 3 and image.shape[0] in (1, 3, 4):
        return ("c", "y", "x")
    raise ValueError(f"Unsupported image shape: {image.shape}")


def add_image(sdata, path, key, coordinate_system, transformation):
    image = read_image(path)
    sdata.images[key] = Image2DModel.parse(
        image,
        dims=image_dims(image),
        transformations={coordinate_system: transformation},
    )
    return image


def add_local_coordinate_system(
    sdata,
    source_coordinate_system,
    local_coordinate_system,
    min_coordinate,
):
    shift_to_local = Translation(
        [-float(min_coordinate[0]), -float(min_coordinate[1])],
        axes=("x", "y"),
    )
    elements = list(sdata.images.values()) + list(sdata.labels.values()) + list(sdata.shapes.values())
    for element in elements:
        original = get_transformation(element, to_coordinate_system=source_coordinate_system)
        set_transformation(
            element,
            Sequence([original, shift_to_local]),
            to_coordinate_system=local_coordinate_system,
        )
    return sdata


def build_block(
    block,
    spaceranger_path,
    fullres_image_path,
    wga_rgb_image_path=None,
    nuclei_polygons_path=None,
    fibres_polygons_path=None,
    background_polygons_path=None,
    roi_metadata_path=None,
    areas_of_interest_path=None,
):
    sdata = visium_hd(
        path=spaceranger_path,
        fullres_image_file=fullres_image_path,
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

    image_key = find_fullres_image_key(sdata, block)
    transformation = get_transformation(sdata.images[image_key], to_coordinate_system=block)

    if wga_rgb_image_path is not None:
        add_image(
            sdata,
            wga_rgb_image_path,
            f"{block}_wga_gfp_dapi_rgb_image",
            block,
            transformation,
        )

    if roi_metadata_path is not None:
        add_geojson_shapes(
            sdata,
            roi_metadata_path,
            "roi_metadata_arvis",
            block,
            transformation,
        )

    if areas_of_interest_path is not None:
        add_geojson_shapes(
            sdata,
            areas_of_interest_path,
            "areas_of_interest_arvis",
            block,
            transformation,
        )

    if nuclei_polygons_path is not None:
        add_geojson_shapes(
            sdata,
            nuclei_polygons_path,
            "nuclei_arvis_shapes",
            block,
            transformation,
        )

    if fibres_polygons_path is not None:
        add_geojson_shapes(
            sdata,
            fibres_polygons_path,
            "fibres_cellpose_shapes",
            block,
            transformation,
        )

    if background_polygons_path is not None:
        add_geojson_shapes(
            sdata,
            background_polygons_path,
            "background_mask_shapes",
            block,
            transformation,
        )

    return sdata


def crop_samples(sdata, block, samples):
    sample_sdatas = {}
    for sample_name, sample_info in samples.items():
        sample_key = sample_info.get("sample_key", f"{block}_{sample_name}")
        x0, y0, x1, y1 = sample_bounds(sample_info)
        sample_sdata = sdata.query.bounding_box(
            axes=["x", "y"],
            min_coordinate=[x0, y0],
            max_coordinate=[x1, y1],
            target_coordinate_system=block,
        )
        add_local_coordinate_system(
            sample_sdata,
            source_coordinate_system=block,
            local_coordinate_system=sample_key,
            min_coordinate=[x0, y0],
        )
        sample_sdata.attrs["block_id"] = block
        sample_sdata.attrs["sample_id"] = sample_name
        sample_sdata.attrs["sample_key"] = sample_key
        sample_sdata.attrs["crop_min_coordinate"] = [x0, y0]
        sample_sdata.attrs["crop_max_coordinate"] = [x1, y1]
        sample_sdata.attrs["rotated_rectangle"] = rotated_rectangle(sample_info)
        sample_sdatas[sample_key] = sample_sdata
    return sample_sdatas


def build_block_and_samples(
    block,
    samples_json,
    spaceranger_path,
    fullres_image_path,
    output_dir,
    wga_rgb_image_path=None,
    nuclei_polygons_path=None,
    fibres_polygons_path=None,
    background_polygons_path=None,
    roi_metadata_path=None,
    areas_of_interest_path=None,
    write_block=True,
    overwrite=True,
):
    samples = load_samples(samples_json, block)

    output_dir = Path(output_dir)
    sdata = build_block(
        block=block,
        spaceranger_path=spaceranger_path,
        fullres_image_path=fullres_image_path,
        wga_rgb_image_path=wga_rgb_image_path,
        nuclei_polygons_path=nuclei_polygons_path,
        fibres_polygons_path=fibres_polygons_path,
        background_polygons_path=background_polygons_path,
        roi_metadata_path=roi_metadata_path,
        areas_of_interest_path=areas_of_interest_path,
    )

    paths = {}
    if write_block:
        paths[block] = write_sdata(sdata, output_dir / f"{block}.zarr", overwrite=overwrite)

    sample_sdatas = crop_samples(sdata, block, samples)
    for sample_key, sample_sdata in sample_sdatas.items():
        paths[sample_key] = write_sdata(
            sample_sdata,
            output_dir / f"{sample_key}.zarr",
            overwrite=overwrite,
        )

    return sdata, sample_sdatas, paths
