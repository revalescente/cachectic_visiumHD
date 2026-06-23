import numpy as np
from cellpose import models


def segment_image(image: np.ndarray, scale: int = 1, gpu: bool = False) -> np.ndarray:
    if scale < 1:
        raise ValueError("scale must be greater than or equal to 1")

    model = models.CellposeModel(gpu=gpu)

    original_shape = image.shape[:2]
    scaled_image = image[::scale, ::scale, ...]

    masks, _, _ = model.eval(
        scaled_image,
        batch_size=32,
        flow_threshold=0.4,
        cellprob_threshold=0.0,
        normalize={"tile_norm_blocksize": 0},
    )

    if scale == 1:
        return masks

    masks = np.repeat(np.repeat(masks, scale, axis=0), scale, axis=1)
    return masks[: original_shape[0], : original_shape[1]]
