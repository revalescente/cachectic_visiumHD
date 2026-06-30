from pathlib import Path


def geojson_path(geojson_dir, block, sample):
    return Path(geojson_dir) / f"{block}_{sample}.geojson"


def translate_coords(coords, dx, dy):
    if isinstance(coords[0], (int, float)):
        return [coords[0] + dx, coords[1] + dy, *coords[2:]]
    return [translate_coords(part, dx, dy) for part in coords]
