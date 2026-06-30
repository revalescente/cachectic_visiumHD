from pathlib import Path
from typing import List, Optional

import typer

app = typer.Typer(
    help=(
        "Cachectic Visium HD pipeline tools: segmentation, polygon conversion, "
        "sample split/merge utilities, and SpatialData construction."
    )
)
utils_app = typer.Typer(
    help=(
        "Utilities for moving between block-level and sample-level files. "
        "TIFF commands split/merge images; GeoJSON commands merge sample annotations."
    )
)
spatial_app = typer.Typer(
    help=(
        "SpatialData commands. Build block/sample zarr stores, merge saved zarr "
        "objects, and run quick visual QC."
    )
)

app.add_typer(utils_app, name="utils")
app.add_typer(spatial_app, name="spatial-data")


def _parse_scales(value: str) -> tuple[int, ...]:
    return tuple(int(scale.strip()) for scale in value.split(",") if scale.strip())


@app.command(
    "background",
    help=(
        "Create background regions from an H&E/RGB TIFF. If OUTPUT ends with "
        ".geojson/.json, writes background polygons; otherwise writes a TIFF mask."
    ),
)
def background(
    image: Path = typer.Option(..., help="Input RGB TIFF image."),
    output: Path = typer.Option(..., help="Output GeoJSON or TIFF mask."),
):
    import tifffile as tiff

    from src.background import background_mask, save_background_geojson

    img = tiff.imread(image)
    mask = background_mask(img)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() in {".geojson", ".json"}:
        save_background_geojson(mask, output)
    else:
        tiff.imwrite(output, (mask.astype("uint8") * 255), compression="zlib")


@app.command(
    "polygons",
    help=(
        "Convert a labeled segmentation TIFF into GeoJSON polygons. Non-zero label "
        "values become polygon features with label_id properties."
    ),
)
def polygons(
    image_mask: Path = typer.Option(..., help="Input labeled TIFF mask."),
    output: Path = typer.Option(..., help="Output GeoJSON polygons."),
):
    import tifffile as tiff

    from src.polygons import save_gdf, to_polygons

    mask = tiff.imread(image_mask)
    save_gdf(to_polygons(mask), output)


@app.command(
    "segmentation",
    help=(
        "Run the multiscale Cellpose segmentation workflow on all images in a "
        "directory. This is the fuller segmentation command with scale merging "
        "and optional rescue tiles."
    ),
)
def segmentation(
    input_dir: Path = typer.Option(Path("INPUTS/H&E"), help="Input image directory."),
    output_dir: Path = typer.Option(Path("OUTPUTS/CellPose"), help="Output mask directory."),
    image_ext: str = typer.Option(".tif", help="Input image extension."),
    scales: str = typer.Option("4,3,2,1", help="Comma-separated multiscale factors."),
    primary_scale: int = typer.Option(4, help="First scale used to initialize the merged mask."),
    rescue_scale: int = typer.Option(2, help="Scale used for rescue tile segmentation."),
    rescue_tiles: bool = typer.Option(True, "--rescue-tiles/--no-rescue-tiles", help="Run extra segmentation on mostly empty tiles."),
    tile_size: int = typer.Option(3072, help="Tile size for rescue segmentation."),
    tile_overlap: int = typer.Option(128, help="Overlap between rescue tiles."),
    min_empty_fraction: float = typer.Option(0.40, help="Minimum empty fraction for a tile to be rescued."),
    min_pixels: int = typer.Option(64, help="Minimum object size retained when merging labels."),
    max_overlap: float = typer.Option(0.20, help="Maximum overlap allowed when adding candidate labels."),
    max_rescue_tiles: int = typer.Option(0, help="Maximum rescue tiles; 0 means no limit."),
    preprocess: bool = typer.Option(True, "--preprocess/--no-preprocess", help="Percentile-normalize image channels before Cellpose."),
    gpu: bool = typer.Option(False, help="Use GPU for Cellpose if available."),
    batch_size: int = typer.Option(32, help="Cellpose batch size."),
    flow_threshold: float = typer.Option(0.4, help="Cellpose flow threshold."),
    cellprob_threshold: float = typer.Option(0.0, help="Cellpose cell probability threshold."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show tqdm progress bars."),
    output_suffix: str = typer.Option("", help="Suffix appended before output mask extension."),
):
    import numpy as np
    from cellpose import io, models
    from natsort import natsorted
    from tqdm import tqdm

    from src.segmentation import segment_image

    image_ext = image_ext if image_ext.startswith(".") else f".{image_ext}"
    files = natsorted(
        [
            path
            for path in input_dir.glob(f"*{image_ext}")
            if "_masks" not in path.name and "_flows" not in path.name
        ]
    )
    if not files:
        raise FileNotFoundError(f"no image files with extension {image_ext} found in {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    masks_ext = ".png" if image_ext == ".png" else ".tif"
    model = models.CellposeModel(gpu=gpu, pretrained_model="cpdino")
    max_tiles = max_rescue_tiles or None

    for path in tqdm(files):
        image = np.array(io.imread(path))
        mask = segment_image(
            image,
            scales=_parse_scales(scales),
            primary_scale=primary_scale,
            rescue_scale=rescue_scale,
            rescue_tiles=rescue_tiles,
            tile_size=tile_size,
            tile_overlap=tile_overlap,
            min_empty_fraction=min_empty_fraction,
            min_pixels=min_pixels,
            max_overlap=max_overlap,
            preprocess=preprocess,
            gpu=gpu,
            model=model,
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            progress=progress,
            max_rescue_tiles=max_tiles,
        )
        io.imsave(output_dir / f"{path.stem}{output_suffix}{masks_ext}", mask)


@app.command(
    "segmentation-simple",
    help=(
        "Run the simple Cellpose wrapper on all images in a directory. Useful for "
        "quick single-scale masks without multiscale merge/rescue logic."
    ),
)
def segmentation_simple(
    input_dir: Path = typer.Option(..., help="Input image directory."),
    output_dir: Path = typer.Option(..., help="Output mask directory."),
    image_ext: str = typer.Option(".tif", help="Input image extension."),
    scale: int = typer.Option(1, help="Downsampling factor before segmentation."),
    gpu: bool = typer.Option(False, help="Use GPU for Cellpose if available."),
):
    import numpy as np
    from cellpose import io
    from natsort import natsorted
    from tqdm import tqdm

    from src.segmentation_simple import segment_image

    image_ext = image_ext if image_ext.startswith(".") else f".{image_ext}"
    files = natsorted(
        [
            path
            for path in input_dir.glob(f"*{image_ext}")
            if "_masks" not in path.name and "_flows" not in path.name
        ]
    )
    if not files:
        raise FileNotFoundError(f"no image files with extension {image_ext} found in {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    masks_ext = ".png" if image_ext == ".png" else ".tif"
    for path in tqdm(files):
        image = np.array(io.imread(path))
        io.imsave(output_dir / f"{path.stem}_masks{masks_ext}", segment_image(image, scale=scale, gpu=gpu))


@utils_app.command(
    "tiff-split",
    help=(
        "Split a block-level TIFF into one TIFF per sample using polygons from "
        "the sample GeoJSON file."
    ),
)
def utils_tiff_split(
    block: str = typer.Argument(..., help="Block name, e.g. blocco1."),
    samples_json: Path = typer.Option(..., "--json", help="Napari-style sample GeoJSON file."),
    image: Path = typer.Option(..., help="Block-level TIFF image to split."),
    output_dir: Path = typer.Option(..., help="Directory for sample TIFF outputs."),
):
    from src.geojson_parser import load_samples
    from src.utils import split_tiff

    split_tiff(block, load_samples(samples_json, block), image, output_dir)


@utils_app.command(
    "tiff-merge",
    help=(
        "Merge sample-level TIFF files back into block coordinates using the "
        "sample GeoJSON polygons."
    ),
)
def utils_tiff_merge(
    block: str = typer.Argument(..., help="Block name, e.g. blocco1."),
    samples_json: Path = typer.Option(..., "--json", help="Napari-style sample GeoJSON file."),
    images_dir: Path = typer.Option(..., help="Directory containing sample TIFF files."),
    output: Path = typer.Option(..., help="Output block-level TIFF."),
):
    from src.geojson_parser import load_samples
    from src.utils import merge_tiff

    merge_tiff(block, load_samples(samples_json, block), images_dir, output)


@utils_app.command(
    "geojson-merge",
    help=(
        "Merge sample-level GeoJSON annotations into block coordinates using the "
        "sample GeoJSON polygons."
    ),
)
def utils_geojson_merge(
    block: str = typer.Argument(..., help="Block name, e.g. blocco1."),
    samples_json: Path = typer.Option(..., "--json", help="Napari-style sample GeoJSON file."),
    geojson_dir: Path = typer.Option(..., help="Directory containing sample-level GeoJSON files."),
    output: Path = typer.Option(..., help="Output block-level GeoJSON."),
):
    from src.geojson_parser import load_samples
    from src.utils import merge_geojson

    merge_geojson(block, load_samples(samples_json, block), geojson_dir, output)


@spatial_app.command(
    "build",
    help=(
        "Build SpatialData zarr stores from SpaceRanger output, H&E full-res image, "
        "optional WGA/GFP/DAPI image, segmentation polygons, background polygons, "
        "and sample GeoJSON polygons. Saves the full block and per-sample crops."
    ),
)
def spatial_build(
    block: str = typer.Argument(..., help="Block name, e.g. blocco1."),
    samples_json: Path = typer.Option(..., help="Napari-style sample GeoJSON file."),
    spaceranger: Path = typer.Option(..., help="SpaceRanger outs directory."),
    fullres_image: Path = typer.Option(..., help="Full-resolution H&E image aligned to SpaceRanger."),
    output_dir: Path = typer.Option(..., help="Directory for SpatialData zarr outputs."),
    wga_rgb_image: Optional[Path] = typer.Option(None, help="Optional WGA/GFP/DAPI RGB image; use H&E placeholder until aligned."),
    nuclei_polygons: Optional[Path] = typer.Option(None, help="GeoJSON polygons from ARVIS nuclei mask."),
    fibres_polygons: Optional[Path] = typer.Option(None, help="GeoJSON polygons from Cellpose fibre mask."),
    background_polygons: Optional[Path] = typer.Option(None, help="GeoJSON background polygons."),
    roi_metadata: Optional[Path] = typer.Option(None, help="Block-level ROI/intissue GeoJSON annotations."),
    areas_of_interest: Optional[Path] = typer.Option(None, help="Optional additional areas-of-interest GeoJSON."),
    write_block: bool = typer.Option(True, "--write-block/--no-block", help="Write the full block zarr in addition to sample zarr stores."),
    overwrite: bool = typer.Option(True, "--overwrite/--no-overwrite", help="Overwrite existing zarr outputs."),
):
    from src.spatial_data.build import build_block_and_samples

    _, _, paths = build_block_and_samples(
        block=block,
        samples_json=samples_json,
        spaceranger_path=spaceranger,
        fullres_image_path=fullres_image,
        wga_rgb_image_path=wga_rgb_image,
        output_dir=output_dir,
        nuclei_polygons_path=nuclei_polygons,
        fibres_polygons_path=fibres_polygons,
        background_polygons_path=background_polygons,
        roi_metadata_path=roi_metadata,
        areas_of_interest_path=areas_of_interest,
        write_block=write_block,
        overwrite=overwrite,
    )
    for key, path in paths.items():
        typer.echo(f"{key} {path}")


@spatial_app.command(
    "merge",
    help=(
        "Merge multiple saved SpatialData zarr stores into one object. Element "
        "keys are prefixed by the provided sample names."
    ),
)
def spatial_merge(
    output: Path = typer.Option(..., help="Output merged SpatialData zarr path."),
    samples: List[str] = typer.Argument(..., help="Pairs like name=path.zarr"),
    overwrite: bool = typer.Option(True, "--overwrite/--no-overwrite"),
):
    from src.spatial_data.merge import merge_and_write

    paths_by_name = {}
    for sample in samples:
        key, path = sample.split("=", maxsplit=1)
        paths_by_name[key] = path
    merge_and_write(paths_by_name, output, overwrite=overwrite)
    typer.echo(output)


@spatial_app.command(
    "qc",
    help=(
        "Summarize and optionally plot a saved SpatialData zarr object. Use this "
        "to inspect image, label, shape, and table keys after build/merge."
    ),
)
def spatial_qc(
    zarr: Optional[Path] = typer.Option(None, help="SpatialData zarr path to inspect."),
    coordinate_system: Optional[str] = typer.Option(None, help="Coordinate system to plot."),
    image_key: Optional[str] = typer.Option(None, help="Image key to render; guessed when omitted."),
    labels_key: Optional[List[str]] = typer.Option(None, help="Label key to render; can be repeated."),
    shapes_key: Optional[List[str]] = typer.Option(None, help="Shape key to render; can be repeated."),
    summary_only: bool = typer.Option(False, help="Print keys/attrs without plotting."),
    output_prefix: Optional[str] = typer.Option(None, help="If provided, save plot as '<prefix>_zarr.png'."),
):
    import matplotlib.pyplot as plt

    from src.spatial_data.qc import plot_zarr, summarize_zarr

    if zarr is None:
        raise typer.BadParameter("--zarr is required")

    summarize_zarr(zarr)
    if summary_only:
        return
    output = f"{output_prefix}_zarr.png" if output_prefix else None
    plot_zarr(
        zarr,
        coordinate_system=coordinate_system,
        image_key=image_key,
        labels_keys=labels_key,
        shapes_keys=shapes_key,
        output=output,
    )
    if output_prefix is None:
        plt.show()


if __name__ == "__main__":
    app()
