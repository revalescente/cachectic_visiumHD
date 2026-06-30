# Cellpose Segmentation

Run multiscale Cellpose segmentation with optional tile rescue on all images in a directory.

## CLI

```bash
python -m src.segmentation \
  --input-dir data/H\&E \
  --output-dir data/CellPose \
  --scales 4,3,2,1 \
  --primary-scale 4 \
  --gpu
```

## Import

```python
from src.segmentation import segment_image, preprocess_image
```
