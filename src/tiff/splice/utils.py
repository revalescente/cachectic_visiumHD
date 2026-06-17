from pathlib import Path


def tile_path(images_dir, block, sample):
    return Path(images_dir) / f"{block}_{sample}.tiff"
