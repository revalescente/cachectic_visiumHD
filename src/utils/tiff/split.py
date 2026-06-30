from pathlib import Path

import tifffile as tiff

from src.geojson_parser import sample_bounds

from .utils import tile_path


def split_tiff(block, samples, image_path, output_dir):
    image = tiff.imread(image_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for sample, info in samples.items():
        x0, y0, x1, y1 = sample_bounds(info)
        tile = image[y0:y1, x0:x1]
        tiff.imwrite(tile_path(output_dir, block, sample), tile, compression="zlib")
