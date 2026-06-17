import argparse
import json

from .merge import merge_tiff
from .split import split_tiff

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["split", "merge"])
    parser.add_argument("block")
    parser.add_argument("--samples-json")
    parser.add_argument("--image")
    parser.add_argument("--images-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--output")
    args = parser.parse_args()

    with open(args.samples_json) as f:
        samples = json.load(f)[args.block]

    if args.mode == "split":
        split_tiff(args.block, samples, args.image, args.output_dir)
    else:
        merge_tiff(args.block, samples, args.images_dir, args.output)


# EXAMPLES
# python split_merge_tiff.py split blocco1 \
#     --image data/full_blocco1.tiff \
#     --output-dir data/split_blocco1
#
# python split_merge_tiff.py merge blocco1 \
#     --images-dir data/split_blocco1 \
#     --output data/merged_blocco1.tiff
