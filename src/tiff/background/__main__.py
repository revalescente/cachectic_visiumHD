import argparse

import tifffile as tiff

from .background import background_mask


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--image")
    parser.add_argument("--output")
    args = parser.parse_args()

    image = tiff.imread(args.image)
    mask = background_mask(image)
    tiff.imwrite(args.output, (mask.astype("uint8") * 255), compression="zlib")


# EXAMPLE
# python -m src.tiff.background \
#     --image data/H&E/Project_BLOCCO_1_ch00.tif \
#     --output data/H&E/Project_BLOCCO_1_background_mask.tif
