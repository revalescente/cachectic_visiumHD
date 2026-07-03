import shutil
from pathlib import Path

from spatialdata import SpatialData, read_zarr


def write_sdata(sdata, path, overwrite=True):
    path = Path(path)
    if path.exists() and overwrite:
        shutil.rmtree(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    sdata.write(path)
    return path


def merge_sdata_objects(paths_by_name):
    merged = SpatialData()
    source_metadata = {}
    for sample_key, path in paths_by_name.items():
        sdata = read_zarr(path)
        source_metadata[sample_key] = {
            "path": str(path),
            "attrs": dict(sdata.attrs),
        }
        for key, element in sdata.images.items():
            merged.images[f"{sample_key}__{key}"] = element
        for key, element in sdata.labels.items():
            merged.labels[f"{sample_key}__{key}"] = element
        for key, element in sdata.shapes.items():
            merged.shapes[f"{sample_key}__{key}"] = element
        for key, table in sdata.tables.items():
            table.obs["source_sample"] = sample_key
            table.obs["source_table"] = key
            if "block_id" in sdata.attrs:
                table.obs["source_block"] = sdata.attrs["block_id"]
            merged.tables[f"{sample_key}__{key}"] = table
    merged.attrs["source_samples"] = list(paths_by_name)
    merged.attrs["source_metadata"] = source_metadata
    return merged


def merge_and_write(paths_by_name, output, overwrite=True):
    merged = merge_sdata_objects(paths_by_name)
    write_sdata(merged, output, overwrite=overwrite)
    return merged
