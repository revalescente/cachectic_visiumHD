import argparse
import json
from pathlib import Path

import cv2
import numpy as np
import tifffile as tiff


def to_gray(image):
    image = np.squeeze(image)
    if image.ndim == 3:
        if image.shape[-1] in (3, 4):
            image = image[..., :3].mean(axis=2)
        elif image.shape[0] in (3, 4):
            image = image[:3].mean(axis=0)
    return image.astype(np.float32)


def to_uint8(image):
    image = to_gray(image)
    low, high = np.percentile(image, (1, 99))
    if high <= low:
        return np.zeros(image.shape, dtype=np.uint8)
    image = np.clip((image - low) / (high - low), 0, 1)
    return (image * 255).astype(np.uint8)


def resize(image, scale):
    return cv2.resize(image, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA)


def resize_to_max(image, max_size):
    h, w = image.shape[:2]
    scale = min(max_size / max(h, w), 1)
    return resize(image, scale), scale


def foreground_mask(image):
    gray = cv2.GaussianBlur(to_uint8(image), (0, 0), 3)
    threshold, _ = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    h, w = gray.shape
    c = max(min(h, w) // 20, 20)
    corners = np.concatenate(
        [
            gray[:c, :c].ravel(),
            gray[:c, -c:].ravel(),
            gray[-c:, :c].ravel(),
            gray[-c:, -c:].ravel(),
        ]
    )

    if np.median(corners) > np.median(gray):
        mask = gray < threshold
    else:
        mask = gray > threshold

    mask = mask.astype(np.uint8) * 255
    kernel = np.ones((9, 9), np.uint8)
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    return mask


def center_of_mask(mask):
    ys, xs = np.nonzero(mask)
    return np.array([xs.mean(), ys.mean()], dtype=np.float32)


def matrix_around_centers(moving_center, fixed_center, scale, angle_degrees):
    angle = np.deg2rad(angle_degrees)
    c = np.cos(angle) * scale
    s = np.sin(angle) * scale
    matrix = np.array([[c, -s, 0], [s, c, 0]], dtype=np.float32)
    matrix[:, 2] = fixed_center - matrix[:, :2] @ moving_center
    return matrix


def warp(image, matrix, shape, interpolation=cv2.INTER_LINEAR):
    h, w = shape
    return cv2.warpAffine(
        image,
        np.asarray(matrix, dtype=np.float32),
        (int(w), int(h)),
        flags=interpolation,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0,
    )


def overlap_score(fixed, moving):
    fixed = fixed > 0
    moving = moving > 0
    union = np.logical_or(fixed, moving).sum()
    if union == 0:
        return 0
    return np.logical_and(fixed, moving).sum() / union


def add_translation(matrix, dx, dy):
    matrix = matrix.copy()
    matrix[0, 2] += dx
    matrix[1, 2] += dy
    return matrix


def translation_by_phase(fixed, moving):
    fixed = fixed.astype(np.float32) / 255
    moving = moving.astype(np.float32) / 255
    window = cv2.createHanningWindow((fixed.shape[1], fixed.shape[0]), cv2.CV_32F)
    shift, response = cv2.phaseCorrelate(fixed * window, moving * window)
    return float(shift[0]), float(shift[1]), float(response)


def refine_candidates(fixed_mask, moving_mask, angles, scales):
    fixed_center = center_of_mask(fixed_mask)
    moving_center = center_of_mask(moving_mask)
    best = None

    for scale in scales:
        for angle in angles:
            matrix = matrix_around_centers(moving_center, fixed_center, scale, angle)
            moved = warp(
                moving_mask,
                matrix,
                fixed_mask.shape,
                interpolation=cv2.INTER_NEAREST,
            )

            dx, dy, response = translation_by_phase(fixed_mask, moved)
            matrix = add_translation(matrix, dx, dy)
            moved = warp(
                moving_mask,
                matrix,
                fixed_mask.shape,
                interpolation=cv2.INTER_NEAREST,
            )
            score = overlap_score(fixed_mask, moved) + response

            if best is None or score > best["score"]:
                best = {
                    "matrix": matrix,
                    "scale": float(scale),
                    "angle": float(angle),
                    "dx": dx,
                    "dy": dy,
                    "response": response,
                    "overlap": float(overlap_score(fixed_mask, moved)),
                    "score": float(score),
                }

    return best


def estimate_similarity(
    fixed_image,
    moving_image,
    max_size=1800,
    min_scale=0.5,
    max_scale=3.0,
    scale_steps=31,
    angle_step=5,
):
    fixed_gray, fixed_resize = resize_to_max(to_gray(fixed_image), max_size)
    moving_gray = resize(to_gray(moving_image), fixed_resize)

    fixed_mask = foreground_mask(fixed_gray)
    moving_mask = foreground_mask(moving_gray)

    angles = np.arange(-180, 180, angle_step)
    scales = np.linspace(min_scale, max_scale, scale_steps)
    coarse = refine_candidates(fixed_mask, moving_mask, angles, scales)

    angles = np.arange(coarse["angle"] - angle_step, coarse["angle"] + angle_step + 0.5, 1)
    scale_delta = (max_scale - min_scale) / max(scale_steps - 1, 1)
    scales = np.linspace(
        max(min_scale, coarse["scale"] - scale_delta),
        min(max_scale, coarse["scale"] + scale_delta),
        11,
    )
    fine = refine_candidates(fixed_mask, moving_mask, angles, scales)

    matrix = fine["matrix"].copy()
    matrix[:, 2] /= fixed_resize
    debug = {
        "fixed_mask": fixed_mask,
        "moving_mask": moving_mask,
        "aligned_moving_mask": warp(
            moving_mask,
            fine["matrix"],
            fixed_mask.shape,
            interpolation=cv2.INTER_NEAREST,
        ),
    }
    return matrix, fine, debug


def save_debug_masks(output_dir, block, debug):
    debug_dir = Path(output_dir) / "debug_masks"
    debug_dir.mkdir(parents=True, exist_ok=True)

    for name, image in debug.items():
        cv2.imwrite(
            str(debug_dir / f"Project_BLOCCO_{block}_{name}.png"),
            image,
        )


def align_block(
    block,
    fixed_dir,
    moving_dir,
    output_dir,
    ref_channel="ch02",
    channels=("ch00", "ch01", "ch02"),
    max_size=1800,
    min_scale=0.5,
    max_scale=3.0,
    scale_steps=31,
    angle_step=5,
):
    fixed_path = Path(fixed_dir) / f"Project_BLOCCO_{block}_ch00.tif"
    moving_ref_path = Path(moving_dir) / f"Project_BLOCCO_{block}_{ref_channel}.tif"

    fixed = tiff.imread(fixed_path)
    moving_ref = tiff.imread(moving_ref_path)

    matrix, info, debug = estimate_similarity(
        fixed,
        moving_ref,
        max_size=max_size,
        min_scale=min_scale,
        max_scale=max_scale,
        scale_steps=scale_steps,
        angle_step=angle_step,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_debug_masks(output_dir, block, debug)

    fixed_shape = to_gray(fixed).shape
    for channel in channels:
        moving_path = Path(moving_dir) / f"Project_BLOCCO_{block}_{channel}.tif"
        moving = tiff.imread(moving_path)
        aligned = warp(moving, matrix, fixed_shape)
        tiff.imwrite(
            output_dir / f"Project_BLOCCO_{block}_{channel}.tif",
            aligned,
            compression="zlib",
        )

    with open(output_dir / f"Project_BLOCCO_{block}_classic_similarity.json", "w") as f:
        info_to_save = {k: v for k, v in info.items() if k != "matrix"}
        json.dump(
            {
                "fixed": str(fixed_path),
                "moving_reference": str(moving_ref_path),
                "matrix_moving_to_fixed": matrix.tolist(),
                "info": info_to_save,
            },
            f,
            indent=2,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("block")
    parser.add_argument("--fixed-dir", default="data/H&E")
    parser.add_argument("--moving-dir", default="data/WGA_GFP_DAPI")
    parser.add_argument("--output-dir", default="data/WGA_GFP_DAPI_aligned_classic")
    parser.add_argument("--ref-channel", default="ch02")
    parser.add_argument("--max-size", type=int, default=1800)
    parser.add_argument("--min-scale", type=float, default=0.5)
    parser.add_argument("--max-scale", type=float, default=3.0)
    parser.add_argument("--scale-steps", type=int, default=31)
    parser.add_argument("--angle-step", type=float, default=5)
    args = parser.parse_args()

    align_block(
        block=args.block,
        fixed_dir=args.fixed_dir,
        moving_dir=args.moving_dir,
        output_dir=args.output_dir,
        ref_channel=args.ref_channel,
        max_size=args.max_size,
        min_scale=args.min_scale,
        max_scale=args.max_scale,
        scale_steps=args.scale_steps,
        angle_step=args.angle_step,
    )


if __name__ == "__main__":
    main()


# EXAMPLES
#
# Classic similarity alignment: rotation + translation + scale
# python -m src.alignment.classic_similarity 1 \
#     --ref-channel ch02 \
#     --output-dir data/WGA_GFP_DAPI_aligned_classic
#
# Wider scale search
# python -m src.alignment.classic_similarity 1 \
#     --min-scale 0.2 \
#     --max-scale 4.0
