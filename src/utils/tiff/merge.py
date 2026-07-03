from pathlib import Path

import numpy as np
import tifffile as tiff

from src.geojson_parser import sample_bounds

from .utils import tile_path


def merge_tiff(block, samples, images_dir, output_path):
    width = max(sample_bounds(info)[2] for info in samples.values())
    height = max(sample_bounds(info)[3] for info in samples.values())

    first_sample = next(iter(samples))
    first_image = tiff.imread(tile_path(images_dir, block, first_sample))
    canvas = np.zeros((height, width, *first_image.shape[2:]), dtype=first_image.dtype)

    for sample, info in samples.items():
        x0, y0, x1, y1 = sample_bounds(info)
        expected_h = y1 - y0
        expected_w = x1 - x0

        image = tiff.imread(tile_path(images_dir, block, sample))[
            :expected_h, :expected_w
        ]
        canvas[y0:y1, x0:x1] = image

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tiff.imwrite(output_path, canvas, compression="zlib")
