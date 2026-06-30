import json
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.features
import shapely.geometry
from spatialdata.models import Labels2DModel
from spatialdata.models import ShapesModel
from spatialdata.transformations import (
    Sequence,
    Translation,
    get_transformation,
    set_transformation,
)
from spatialdata_io import visium_hd


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


def read_label_raster(path):
    with rasterio.open(path) as src:
        return src.read(1)


def mask_to_polygons(mask):
    geometries = []
    labels = []
    mask = np.asarray(mask)
    for geometry, label_id in rasterio.features.shapes(
        mask.astype(np.int32),
        mask=mask > 0,
        connectivity=8,
    ):
        geometries.append(shapely.geometry.shape(geometry))
        labels.append(int(label_id))

    gdf = gpd.GeoDataFrame({"label_id": labels, "geometry": geometries}, index=labels)
    if len(labels) != len(set(labels)):
        gdf = gdf.dissolve(by="label_id")
        gdf["label_id"] = gdf.index.astype(int)
    return gdf.set_crs(None, allow_override=True)


def add_geojson_shapes(sdata, path, key, coordinate_system, transformation):
    gdf = gpd.read_file(path).set_crs(None, allow_override=True)
    gdf.attrs = {}
    sdata.shapes[key] = ShapesModel.parse(
        gdf,
        transformations={coordinate_system: transformation},
    )
    return gdf


def add_label_mask(
    sdata,
    path,
    labels_key,
    shapes_key,
    coordinate_system,
    transformation,
    make_polygons=True,
):
    mask = read_label_raster(path)
    sdata.labels[labels_key] = Labels2DModel.parse(
        mask,
        dims=("y", "x"),
        transformations={coordinate_system: transformation},
    )
    if make_polygons:
        polygons = mask_to_polygons(mask)
        sdata.shapes[shapes_key] = ShapesModel.parse(
            polygons,
            transformations={coordinate_system: transformation},
        )
    return mask


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
    nuclei_labels_path,
    fibres_labels_path,
    roi_metadata_path=None,
    areas_of_interest_path=None,
    make_fibre_polygons=True,
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

    nuclei_mask = read_label_raster(nuclei_labels_path)
    sdata.labels["nuclei_arvis_labels"] = Labels2DModel.parse(
        nuclei_mask,
        dims=("y", "x"),
        transformations={block: transformation},
    )

    add_label_mask(
        sdata,
        fibres_labels_path,
        "fibres_cellpose_labels",
        "fibres_cellpose_shapes",
        block,
        transformation,
        make_polygons=make_fibre_polygons,
    )

    return sdata


def crop_samples(sdata, block, samples):
    sample_sdatas = {}
    for sample_name, sample_info in samples.items():
        sample_key = sample_info.get("sample_key", f"{block}_{sample_name}")
        sample_sdata = sdata.query.bounding_box(
            axes=["x", "y"],
            min_coordinate=sample_info["min_coordinate"],
            max_coordinate=sample_info["max_coordinate"],
            target_coordinate_system=block,
        )
        add_local_coordinate_system(
            sample_sdata,
            source_coordinate_system=block,
            local_coordinate_system=sample_key,
            min_coordinate=sample_info["min_coordinate"],
        )
        sample_sdata.attrs["block_id"] = block
        sample_sdata.attrs["sample_id"] = sample_name
        sample_sdata.attrs["sample_key"] = sample_key
        sample_sdata.attrs["crop_min_coordinate"] = sample_info["min_coordinate"]
        sample_sdata.attrs["crop_max_coordinate"] = sample_info["max_coordinate"]
        sample_sdatas[sample_key] = sample_sdata
    return sample_sdatas


def build_block_and_samples(
    block,
    samples_json,
    spaceranger_path,
    fullres_image_path,
    nuclei_labels_path,
    fibres_labels_path,
    output_dir,
    roi_metadata_path=None,
    areas_of_interest_path=None,
    make_fibre_polygons=True,
    write_block=True,
    overwrite=True,
):
    with open(samples_json) as f:
        samples = json.load(f)[block]

    output_dir = Path(output_dir)
    sdata = build_block(
        block=block,
        spaceranger_path=spaceranger_path,
        fullres_image_path=fullres_image_path,
        nuclei_labels_path=nuclei_labels_path,
        fibres_labels_path=fibres_labels_path,
        roi_metadata_path=roi_metadata_path,
        areas_of_interest_path=areas_of_interest_path,
        make_fibre_polygons=make_fibre_polygons,
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
