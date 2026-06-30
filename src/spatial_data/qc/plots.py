import json

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.patches import Rectangle
from rasterio.windows import Window


def image_size(path):
    with rasterio.open(path) as src:
        return src.width, src.height


def read_rgb_overview(path, max_size=2200):
    with rasterio.open(path) as src:
        scale = min(max_size / src.width, max_size / src.height, 1.0)
        out_w = max(1, int(src.width * scale))
        out_h = max(1, int(src.height * scale))
        indexes = [1, 2, 3] if src.count >= 3 else [1]
        arr = src.read(indexes, out_shape=(len(indexes), out_h, out_w))
        arr = np.moveaxis(arr, 0, -1)
    if arr.shape[-1] == 1:
        arr = arr[..., 0]
    return arr


def read_rgb_window(path, x0, y0, width, height):
    with rasterio.open(path) as src:
        x0 = int(max(0, min(x0, src.width - 1)))
        y0 = int(max(0, min(y0, src.height - 1)))
        width = int(max(1, min(width, src.width - x0)))
        height = int(max(1, min(height, src.height - y0)))
        indexes = [1, 2, 3] if src.count >= 3 else [1]
        arr = src.read(indexes, window=Window(x0, y0, width, height))
        arr = np.moveaxis(arr, 0, -1)
    if arr.shape[-1] == 1:
        arr = arr[..., 0]
    return arr, (x0, y0, width, height)


def read_label_window(path, x0, y0, width, height):
    with rasterio.open(path) as src:
        x0 = int(max(0, min(x0, src.width - 1)))
        y0 = int(max(0, min(y0, src.height - 1)))
        width = int(max(1, min(width, src.width - x0)))
        height = int(max(1, min(height, src.height - y0)))
        arr = src.read(1, window=Window(x0, y0, width, height))
    return arr


def sample_bbox(info):
    x0, y0 = info["min_coordinate"]
    x1, y1 = info["max_coordinate"]
    return int(x0), int(y0), int(x1 - x0), int(y1 - y0)


def centered_window(info, size=1400):
    x, y, w, h = sample_bbox(info)
    win_w = min(size, w)
    win_h = min(size, h)
    x0 = x + (w - win_w) // 2
    y0 = y + (h - win_h) // 2
    return int(x0), int(y0), int(win_w), int(win_h)


def plot_block_overview(
    he_path,
    samples,
    roi_metadata_path=None,
    areas_of_interest_path=None,
    max_size=2400,
    ax=None,
):
    overview = read_rgb_overview(he_path, max_size=max_size)
    full_w, full_h = image_size(he_path)
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 13))
    ax.imshow(overview, extent=[0, full_w, full_h, 0])

    if roi_metadata_path is not None:
        roi_metadata = gpd.read_file(roi_metadata_path)
        roi_metadata.boundary.plot(ax=ax, color="lime", linewidth=1.0, alpha=0.9)

    if areas_of_interest_path is not None:
        areas = gpd.read_file(areas_of_interest_path)
        areas.geometry.boundary.plot(ax=ax, color="white", linewidth=0.15, alpha=0.25)

    for sample_name, info in samples.items():
        x, y, w, h = sample_bbox(info)
        ax.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor="red", linewidth=2))
        ax.text(
            x + 120,
            y + 350,
            info.get("sample_key", sample_name),
            color="red",
            fontsize=11,
            weight="bold",
        )

    ax.set_title("Block overview")
    ax.set_xlim(0, full_w)
    ax.set_ylim(full_h, 0)
    ax.set_aspect("equal")
    return ax


def plot_segmentation_window(
    he_path,
    nuclei_labels_path,
    fibres_labels_path,
    x0,
    y0,
    width,
    height,
    title=None,
    ax=None,
    global_coordinates=False,
):
    he, (x0, y0, width, height) = read_rgb_window(he_path, x0, y0, width, height)
    nuclei = read_label_window(nuclei_labels_path, x0, y0, width, height)
    fibres = read_label_window(fibres_labels_path, x0, y0, width, height)
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 8))

    if global_coordinates:
        extent = [x0, x0 + width, y0 + height, y0]
        yy, xx = np.mgrid[y0 : y0 + height, x0 : x0 + width]
    else:
        extent = [0, width, height, 0]
        yy, xx = np.mgrid[0:height, 0:width]

    ax.imshow(he, extent=extent)
    if np.any(nuclei > 0):
        ax.contour(xx, yy, nuclei > 0, levels=[0.5], colors=["cyan"], linewidths=0.5)
    if np.any(fibres > 0):
        ax.contour(xx, yy, fibres > 0, levels=[0.5], colors=["yellow"], linewidths=0.7)

    ax.set_title(title or f"x={x0}, y={y0}, w={width}, h={height}")
    if global_coordinates:
        ax.set_xlim(x0, x0 + width)
        ax.set_ylim(y0 + height, y0)
    else:
        ax.set_xlim(0, width)
        ax.set_ylim(height, 0)
    ax.set_aspect("equal")
    return ax


def plot_sample_windows(
    he_path,
    nuclei_labels_path,
    fibres_labels_path,
    samples,
    size=1400,
):
    fig, axes = plt.subplots(1, len(samples), figsize=(6 * len(samples), 6), constrained_layout=True)
    if len(samples) == 1:
        axes = [axes]
    for ax, (sample_name, info) in zip(axes, samples.items()):
        x0, y0, w, h = centered_window(info, size=size)
        plot_segmentation_window(
            he_path,
            nuclei_labels_path,
            fibres_labels_path,
            x0,
            y0,
            w,
            h,
            title=info.get("sample_key", sample_name),
            ax=ax,
        )
    return fig, axes


def read_samples(path, block):
    with open(path) as f:
        return json.load(f)[block]
