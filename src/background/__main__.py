import argparse
from pathlib import Path

import tifffile as tiff

from .background import background_mask, save_background_geojson


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--image")
    parser.add_argument("--output")
    args = parser.parse_args()

    image = tiff.imread(args.image)
    mask = background_mask(image)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix.lower() in {".geojson", ".json"}:
        save_background_geojson(mask, output)
    else:
        tiff.imwrite(output, (mask.astype("uint8") * 255), compression="zlib")


# EXAMPLE
# python -m src.background \
#     --image data/H&E/Project_BLOCCO_1_ch00.tif \
#     --output data/H&E/Project_BLOCCO_1_background_mask.geojson
