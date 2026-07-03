import argparse
from pathlib import Path

import numpy as np
from cellpose import io, models
from natsort import natsorted
from tqdm import tqdm

from .cellpose import segment_image


def _parse_scales(value: str) -> tuple[int, ...]:
    return tuple(int(scale.strip()) for scale in value.split(",") if scale.strip())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="data/H&E")
    parser.add_argument("--output-dir", default="data/CellPose")
    parser.add_argument("--image-ext", default=".tif")
    parser.add_argument("--scales", default="4,3,2,1")
    parser.add_argument("--primary-scale", type=int, default=4)
    parser.add_argument("--rescue-scale", type=int, default=2)
    parser.add_argument("--no-rescue-tiles", action="store_true")
    parser.add_argument("--tile-size", type=int, default=3072)
    parser.add_argument("--tile-overlap", type=int, default=128)
    parser.add_argument("--min-empty-fraction", type=float, default=0.40)
    parser.add_argument("--min-pixels", type=int, default=64)
    parser.add_argument("--max-overlap", type=float, default=0.20)
    parser.add_argument("--max-rescue-tiles", type=int, default=0)
    parser.add_argument("--no-preprocess", action="store_true")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--flow-threshold", type=float, default=0.4)
    parser.add_argument("--cellprob-threshold", type=float, default=0.0)
    parser.add_argument("--no-progress", action="store_true")
    parser.add_argument("--output-suffix", default="")

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    image_ext = (
        args.image_ext if args.image_ext.startswith(".") else f".{args.image_ext}"
    )

    if not input_dir.exists():
        raise FileNotFoundError(f"input directory does not exist: {input_dir}")

    files = natsorted(
        [
            path
            for path in input_dir.glob(f"*{image_ext}")
            if "_masks" not in path.name and "_flows" not in path.name
        ]
    )
    if len(files) == 0:
        raise FileNotFoundError(
            f"no image files with extension {image_ext} found in {input_dir}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    masks_ext = ".png" if image_ext == ".png" else ".tif"
    scales = _parse_scales(args.scales)
    model = models.CellposeModel(gpu=args.gpu, pretrained_model="cpdino")
    max_rescue_tiles = args.max_rescue_tiles or None

    for path in tqdm(files):
        image = np.array(io.imread(path))
        segmentation = segment_image(
            image,
            scales=scales,
            primary_scale=args.primary_scale,
            rescue_scale=args.rescue_scale,
            rescue_tiles=not args.no_rescue_tiles,
            tile_size=args.tile_size,
            tile_overlap=args.tile_overlap,
            min_empty_fraction=args.min_empty_fraction,
            min_pixels=args.min_pixels,
            max_overlap=args.max_overlap,
            preprocess=not args.no_preprocess,
            gpu=args.gpu,
            model=model,
            batch_size=args.batch_size,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
            progress=not args.no_progress,
            max_rescue_tiles=max_rescue_tiles,
        )
        io.imsave(
            output_dir / f"{path.stem}{args.output_suffix}{masks_ext}",
            segmentation,
        )


# EXAMPLE
# python -m src.segmentation \
#     --input-dir data/H&E \
#     --output-dir data/CellPose \
#     --scales 4,3,2,1 \
#     --primary-scale 4 \
#     --gpu
