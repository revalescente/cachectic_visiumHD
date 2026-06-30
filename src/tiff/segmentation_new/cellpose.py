from __future__ import annotations

import numpy as np
from cellpose import models
from tqdm import tqdm


def preprocess_image(image: np.ndarray, low: float = 1, high: float = 99) -> np.ndarray:
    image = image.astype(np.float32, copy=False)
    out = np.empty_like(image, dtype=np.float32)

    if image.ndim == 2:
        p_low, p_high = np.percentile(image, [low, high])
        out = (image - p_low) / max(p_high - p_low, 1)
        return np.clip(out, 0, 1)

    for channel in range(image.shape[-1]):
        plane = image[..., channel]
        p_low, p_high = np.percentile(plane, [low, high])
        out[..., channel] = (plane - p_low) / max(p_high - p_low, 1)

    return np.clip(out, 0, 1)


def _upsample_labels(mask: np.ndarray, scale: int, shape: tuple[int, int]) -> np.ndarray:
    if scale == 1:
        return mask[: shape[0], : shape[1]]

    mask = np.repeat(np.repeat(mask, scale, axis=0), scale, axis=1)
    return mask[: shape[0], : shape[1]]


def _segment_with_model(
    model: models.CellposeModel,
    image: np.ndarray,
    scale: int,
    batch_size: int,
    flow_threshold: float,
    cellprob_threshold: float,
) -> np.ndarray:
    original_shape = image.shape[:2]
    scaled_image = image[::scale, ::scale, ...]

    masks, _, _ = model.eval(
        scaled_image,
        batch_size=batch_size,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        normalize=True,
    )

    return _upsample_labels(masks, scale, original_shape)


def _merge_labels(
    base: np.ndarray,
    candidate: np.ndarray,
    min_pixels: int,
    max_overlap: float,
    start_id: int | None = None,
) -> np.ndarray:
    out = base.copy()
    next_id = int(out.max()) + 1 if start_id is None else start_id
    occupied = out > 0
    candidate = candidate.astype(np.int64, copy=False)

    sizes = np.bincount(candidate.ravel())
    if len(sizes) <= 1:
        return out

    overlap_sizes = np.bincount(candidate[occupied].ravel(), minlength=len(sizes))
    labels = np.flatnonzero(sizes >= min_pixels)
    labels = labels[labels != 0]

    for label in labels:
        pixels = int(sizes[label])
        overlap = int(overlap_sizes[label]) / pixels

        if overlap > max_overlap:
            continue

        region = candidate == label
        region = region & ~occupied
        if int(region.sum()) < min_pixels:
            continue

        out[region] = next_id
        occupied[region] = True
        next_id += 1

    return out


def _iter_tiles(
    height: int,
    width: int,
    tile_size: int,
    overlap: int,
):
    step = max(tile_size - overlap, 1)
    ys = list(range(0, height, step))
    xs = list(range(0, width, step))
    seen = set()

    for y in ys:
        y0 = min(y, max(height - tile_size, 0))
        y1 = min(y0 + tile_size, height)
        for x in xs:
            x0 = min(x, max(width - tile_size, 0))
            x1 = min(x0 + tile_size, width)
            key = (y0, y1, x0, x1)
            if key in seen:
                continue
            seen.add(key)
            yield y0, y1, x0, x1


def _rescue_unsegmented_tiles(
    model: models.CellposeModel,
    image: np.ndarray,
    mask: np.ndarray,
    scale: int,
    tile_size: int,
    overlap: int,
    min_empty_fraction: float,
    min_pixels: int,
    max_overlap: float,
    batch_size: int,
    flow_threshold: float,
    cellprob_threshold: float,
    progress: bool,
    max_tiles: int | None,
) -> np.ndarray:
    out = mask.copy()
    height, width = out.shape[:2]
    tiles = []
    for y0, y1, x0, x1 in _iter_tiles(height, width, tile_size, overlap):
        tile_mask = out[y0:y1, x0:x1]
        empty_fraction = float((tile_mask == 0).mean())
        if empty_fraction >= min_empty_fraction:
            tiles.append((empty_fraction, y0, y1, x0, x1))

    tiles.sort(reverse=True)
    if max_tiles is not None:
        tiles = tiles[:max_tiles]

    iterator = tqdm(tiles, desc="rescue tiles", leave=False) if progress else tiles
    for _, y0, y1, x0, x1 in iterator:
        tile_mask = out[y0:y1, x0:x1]
        tile_image = image[y0:y1, x0:x1, ...]
        candidate = _segment_with_model(
            model,
            tile_image,
            scale=scale,
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )

        merged_tile = _merge_labels(
            tile_mask,
            candidate,
            min_pixels=min_pixels,
            max_overlap=max_overlap,
            start_id=int(out.max()) + 1,
        )
        out[y0:y1, x0:x1] = merged_tile

    return out


def segment_image(
    image: np.ndarray,
    scales: tuple[int, ...] = (4, 3, 2, 1),
    primary_scale: int = 4,
    rescue_scale: int = 2,
    rescue_tiles: bool = True,
    tile_size: int = 3072,
    tile_overlap: int = 128,
    min_empty_fraction: float = 0.40,
    min_pixels: int = 64,
    max_overlap: float = 0.20,
    preprocess: bool = True,
    gpu: bool = False,
    model: models.CellposeModel | None = None,
    batch_size: int = 32,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    progress: bool = True,
    max_rescue_tiles: int | None = None,
) -> np.ndarray:
    if primary_scale not in scales:
        scales = (primary_scale, *scales)

    scales = tuple(dict.fromkeys(scales))
    if any(scale < 1 for scale in scales):
        raise ValueError("scales must be greater than or equal to 1")

    image_for_model = preprocess_image(image) if preprocess else image
    if model is None:
        model = models.CellposeModel(gpu=gpu, pretrained_model="cpdino")

    ordered_scales = (primary_scale, *[scale for scale in scales if scale != primary_scale])
    merged = None

    scale_iterator = (
        tqdm(ordered_scales, desc="scales", leave=False) if progress else ordered_scales
    )
    for scale in scale_iterator:
        candidate = _segment_with_model(
            model,
            image_for_model,
            scale=scale,
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
        )

        if merged is None:
            merged = candidate.astype(np.uint32, copy=False)
        else:
            merged = _merge_labels(
                merged,
                candidate,
                min_pixels=min_pixels,
                max_overlap=max_overlap,
            )

    if rescue_tiles:
        merged = _rescue_unsegmented_tiles(
            model,
            image_for_model,
            merged,
            scale=rescue_scale,
            tile_size=tile_size,
            overlap=tile_overlap,
            min_empty_fraction=min_empty_fraction,
            min_pixels=min_pixels,
            max_overlap=max_overlap,
            batch_size=batch_size,
            flow_threshold=flow_threshold,
            cellprob_threshold=cellprob_threshold,
            progress=progress,
            max_tiles=max_rescue_tiles,
        )

    return merged.astype(np.uint32, copy=False)
