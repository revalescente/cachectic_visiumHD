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


# EXAMPLES
#
# ARVIS mask -> GeoJSON polygons
# python -m src.polygons \
#     --image-mask data/ARVIS/Project_Blocco_1_ch00.tif \
#     --output data/ARVIS/Project_Blocco_1_ch00.geojson
#
# Cellpose mask -> GeoJSON polygons
# python -m src.polygons \
#     --image-mask data/cellpose/blocco1_cellpose_masks.tiff \
#     --output data/cellpose/blocco1_cellpose_polygons.geojson
