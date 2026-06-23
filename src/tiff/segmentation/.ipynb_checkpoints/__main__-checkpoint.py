import argparse
from pathlib import Path

import numpy as np
from cellpose import io
from natsort import natsorted
from tqdm import tqdm

from .cellpose import segment_image

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--image-ext", default=".tif")
    parser.add_argument("--scale", type=int, default=1)
    parser.add_argument("--gpu", action="store_true")

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

    for path in tqdm(files):
        image = np.array(io.imread(path))
        segmentation = segment_image(image, scale=args.scale, gpu=args.gpu)
        io.imsave(output_dir / f"{path.stem}_masks{masks_ext}", segmentation)


# EXAMPLE
# python -m src.tiff.segmentation \
#     --input-dir data/H&E/Project_BLOCCO_1_ch00 \
#     --output-dir data/cellpose \
#     --image-ext .tif \
#     --scale 4
