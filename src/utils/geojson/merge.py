import copy
import json
from pathlib import Path

from src.geojson_parser import sample_bounds

from .utils import geojson_path, translate_coords


def merge_geojson(block, samples, geojson_dir, output_path):
    features = []

    for sample, info in samples.items():
        x0, y0, _, _ = sample_bounds(info)

        with open(geojson_path(geojson_dir, block, sample)) as f:
            data = json.load(f)

        for feature in data["features"]:
            shifted = copy.deepcopy(feature)
            geometry = shifted.get("geometry")

            if geometry and geometry.get("coordinates"):
                geometry["coordinates"] = translate_coords(
                    geometry["coordinates"],
                    x0,
                    y0,
                )

            shifted.setdefault("properties", {})["sample"] = sample
            shifted.setdefault("properties", {})["block"] = block
            features.append(shifted)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)
