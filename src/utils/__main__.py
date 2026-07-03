import argparse

from src.geojson_parser import load_samples
from .geojson import merge as merge_geojson
from .tiff import merge as merge_tiff
from .tiff import split as split_tiff

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=["tiff", "geojson"])
    parser.add_argument("mode", choices=["split", "merge"])
    parser.add_argument("block")
    parser.add_argument("--json", required=True)
    parser.add_argument("--image")
    parser.add_argument("--images-dir")
    parser.add_argument("--geojson-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--output")
    args = parser.parse_args()

    samples = load_samples(args.json, args.block)

    if args.kind == "tiff" and args.mode == "split":
        split_tiff(args.block, samples, args.image, args.output_dir)
    elif args.kind == "tiff" and args.mode == "merge":
        merge_tiff(args.block, samples, args.images_dir, args.output)
    elif args.kind == "geojson" and args.mode == "merge":
        merge_geojson(args.block, samples, args.geojson_dir, args.output)
    else:
        raise ValueError("geojson split is not implemented")


# EXAMPLES
# python -m src.utils tiff split blocco1 \
#     --json data/samples_he.json \
#     --image data/H\&E/Project_BLOCCO_1_ch00.tif \
#     --output-dir data/H\&E/Project_BLOCCO_1_ch00
#
# python -m src.utils tiff merge blocco1 \
#     --json data/samples_he.json \
#     --images-dir data/ARVIS/arivis_blocco1 \
#     --output data/ARVIS/Project_BLOCCO_1_ch00.tif
#
# python -m src.utils geojson merge blocco1 \
#     --json data/samples_he.json \
#     --geojson-dir data/ARVIS/intissue_GFP_blocco1 \
#     --output data/ARVIS/Project_BLOCCO_1_ch00.geojson
