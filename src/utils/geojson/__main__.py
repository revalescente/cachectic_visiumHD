import argparse

from src.geojson_parser import load_samples
from .merge import merge_geojson

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("block")
    parser.add_argument("--json")
    parser.add_argument("--geojson-dir")
    parser.add_argument("--output")
    args = parser.parse_args()

    samples = load_samples(args.json, args.block)

    merge_geojson(args.block, samples, args.geojson_dir, args.output)


# EXAMPLE
# python -m src.utils geojson merge blocco1 \
#     --json data/samples.json \
#     --geojson-dir data/ARVIS/intissue_GFP_blocco1 \
#     --output data/ARVIS/Project_BLOCCO_1_ch00.geojson
