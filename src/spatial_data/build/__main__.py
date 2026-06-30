import argparse

from .builder import build_block_and_samples

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("block")
    parser.add_argument("--samples-json", required=True)
    parser.add_argument("--spaceranger", required=True)
    parser.add_argument("--fullres-image", required=True)
    parser.add_argument("--nuclei-labels", required=True)
    parser.add_argument("--fibres-labels", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--roi-metadata")
    parser.add_argument("--areas-of-interest")
    parser.add_argument("--no-fibre-polygons", action="store_true")
    parser.add_argument("--no-block", action="store_true")
    parser.add_argument("--no-overwrite", action="store_true")
    args = parser.parse_args()

    _, _, paths = build_block_and_samples(
        block=args.block,
        samples_json=args.samples_json,
        spaceranger_path=args.spaceranger,
        fullres_image_path=args.fullres_image,
        nuclei_labels_path=args.nuclei_labels,
        fibres_labels_path=args.fibres_labels,
        output_dir=args.output_dir,
        roi_metadata_path=args.roi_metadata,
        areas_of_interest_path=args.areas_of_interest,
        make_fibre_polygons=not args.no_fibre_polygons,
        write_block=not args.no_block,
        overwrite=not args.no_overwrite,
    )

    for key, path in paths.items():
        print(key, path)


# EXAMPLE
# python -m src.spatial_data.build blocco1 \
#     --samples-json data/samples_he.json \
#     --spaceranger data/SpaceRanger/blocco1/outs \
#     --fullres-image data/H\&E/Project_BLOCCO_1_ch00.tif \
#     --nuclei-labels data/ARVIS/Project_BLOCCO_1_ch00.tif \
#     --fibres-labels data/CellPose/Project_BLOCCO_1_ch00.tif \
#     --roi-metadata data/ARVIS/Project_BLOCCO_1_ch00.geojson \
#     --output-dir data/spatialdata
