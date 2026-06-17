import argparse

import tifffile as tiff

from .polygons import save_gdf, to_polygons

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--image-mask")
    parser.add_argument("--output")
    args = parser.parse_args()

    image = tiff.imread(args.image_mask)
    polygons = to_polygons(image)
    save_gdf(polygons, args.output)
