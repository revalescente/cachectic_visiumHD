import argparse

from .merge import merge_and_write

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--no-overwrite", action="store_true")
    parser.add_argument("samples", nargs="+")
    args = parser.parse_args()

    paths_by_name = {}
    for sample in args.samples:
        key, path = sample.split("=", maxsplit=1)
        paths_by_name[key] = path

    merge_and_write(
        paths_by_name,
        args.output,
        overwrite=not args.no_overwrite,
    )
    print(args.output)


# EXAMPLE
# python -m src.spatial_data.merge \
#     --output data/spatialdata/samples_merged.zarr \
#     blocco1_c26STAT3=data/spatialdata/blocco1_c26STAT3.zarr \
#     blocco1_sham=data/spatialdata/blocco1_sham.zarr \
#     blocco1_c26foxO=data/spatialdata/blocco1_c26foxO.zarr \
#     blocco2_c26=data/spatialdata/blocco2_c26.zarr
