import argparse

import matplotlib.pyplot as plt

from .plots import plot_block_overview, plot_sample_windows, read_samples
from .zarr import plot_zarr, summarize_zarr

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("block", nargs="?")
    parser.add_argument("--zarr")
    parser.add_argument("--coordinate-system")
    parser.add_argument("--image-key")
    parser.add_argument("--labels-key", action="append")
    parser.add_argument("--shapes-key", action="append")
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--samples-json")
    parser.add_argument("--he")
    parser.add_argument("--nuclei-labels")
    parser.add_argument("--fibres-labels")
    parser.add_argument("--areas-of-interest")
    parser.add_argument("--roi-metadata")
    parser.add_argument("--output-prefix")
    args = parser.parse_args()

    if args.zarr:
        summarize_zarr(args.zarr)
        if not args.summary_only:
            output = f"{args.output_prefix}_zarr.png" if args.output_prefix else None
            plot_zarr(
                args.zarr,
                coordinate_system=args.coordinate_system,
                image_key=args.image_key,
                labels_keys=args.labels_key,
                shapes_keys=args.shapes_key,
                output=output,
            )
            if not args.output_prefix:
                plt.show()
        raise SystemExit(0)

    if not all([args.block, args.samples_json, args.he, args.nuclei_labels, args.fibres_labels]):
        parser.error("raw QC richiede block, --samples-json, --he, --nuclei-labels e --fibres-labels; oppure usa --zarr")

    samples = read_samples(args.samples_json, args.block)

    plot_block_overview(
        args.he,
        samples,
        roi_metadata_path=args.roi_metadata,
        areas_of_interest_path=args.areas_of_interest,
    )
    if args.output_prefix:
        plt.savefig(f"{args.output_prefix}_overview.png", dpi=160)
    else:
        plt.show()

    plot_sample_windows(
        args.he,
        args.nuclei_labels,
        args.fibres_labels,
        samples,
    )
    if args.output_prefix:
        plt.savefig(f"{args.output_prefix}_samples.png", dpi=160)
    else:
        plt.show()


# EXAMPLE
# python -m src.spatial_data.qc --zarr data/spatialdata/blocco1_sham.zarr
#
# python -m src.spatial_data.qc --zarr data/spatialdata/blocco1_sham.zarr \
#     --coordinate-system blocco1_sham \
#     --output-prefix data/spatialdata/qc/blocco1_sham
#
# python -m src.spatial_data.qc blocco1 \
#     --samples-json data/samples_he.json \
#     --he data/H\&E/Project_BLOCCO_1_ch00.tif \
#     --nuclei-labels data/ARVIS/Project_BLOCCO_1_ch00.tif \
#     --fibres-labels data/CellPose/Project_BLOCCO_1_ch00.tif \
#     --areas-of-interest data/ARVIS/Project_Blocco_1_ch00_poli.geojson \
#     --roi-metadata data/ARVIS/Project_BLOCCO_1_ch00.geojson \
#     --output-prefix data/spatialdata/qc/blocco1
