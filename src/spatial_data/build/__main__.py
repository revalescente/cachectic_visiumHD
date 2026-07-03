import argparse

from .builder import build_block_and_samples

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("block")
    parser.add_argument("--samples-json", required=True)
    parser.add_argument("--spaceranger", required=True)
    parser.add_argument("--fullres-image", required=True)
    parser.add_argument("--wga-rgb-image")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--nuclei-polygons")
    parser.add_argument("--fibres-polygons")
    parser.add_argument("--background-polygons")
    parser.add_argument("--roi-metadata")
    parser.add_argument("--areas-of-interest")
    parser.add_argument("--no-block", action="store_true")
    parser.add_argument("--no-overwrite", action="store_true")
    args = parser.parse_args()

    _, _, paths = build_block_and_samples(
        block=args.block,
        samples_json=args.samples_json,
        spaceranger_path=args.spaceranger,
        fullres_image_path=args.fullres_image,
        wga_rgb_image_path=args.wga_rgb_image,
        output_dir=args.output_dir,
        nuclei_polygons_path=args.nuclei_polygons,
        fibres_polygons_path=args.fibres_polygons,
        background_polygons_path=args.background_polygons,
        roi_metadata_path=args.roi_metadata,
        areas_of_interest_path=args.areas_of_interest,
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
#     --wga-rgb-image data/H\&E/Project_BLOCCO_1_ch00.tif \
#     --nuclei-polygons data/ARVIS/Project_BLOCCO_1_ch00_polygons.geojson \
#     --fibres-polygons data/CellPose/Project_BLOCCO_1_ch00_polygons.geojson \
#     --background-polygons data/H\&E/Project_BLOCCO_1_background_mask.geojson \
#     --roi-metadata data/ARVIS/Project_BLOCCO_1_ch00.geojson \
#     --output-dir data/spatialdata
