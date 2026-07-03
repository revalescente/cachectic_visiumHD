import argparse

from src.geojson_parser import load_samples
from .merge import merge_tiff
from .split import split_tiff

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["split", "merge"])
    parser.add_argument("block")
    parser.add_argument("--json")
    parser.add_argument("--image")
    parser.add_argument("--images-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--output")
    args = parser.parse_args()

    samples = load_samples(args.json, args.block)

    if args.mode == "split":
        split_tiff(args.block, samples, args.image, args.output_dir)
    else:
        merge_tiff(args.block, samples, args.images_dir, args.output)


# EXAMPLES
# python -m src.utils tiff split blocco1 \
#    --json data/samples.json \
#     --image "data/H&E/Project_BLOCCO_1_ch00.tif" \
#     --output-dir "data/H&E/Project_BLOCCO_1_ch00"
#
# python -m src.utils tiff merge blocco1 \
#    --json data/samples.json \
#     --images-dir "data/H&E/Project_BLOCCO_1_ch00" \
#     --output "data/H&E/Project_BLOCCO_1_ch00_rec.tif"

# python -m src.utils tiff merge blocco1 \
#    --json data/samples.json \
#     --images-dir "data/ARVIS/arivis_blocco1" \
#     --output "data/ARVIS/Project_BLOCCO_1_ch00.tif"
