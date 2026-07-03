# Cellpose Segmentation Simple

Run a simple Cellpose segmentation on all images in a directory.

## CLI

```bash
python -m src.segmentation_simple \
  --input-dir data/H\&E \
  --output-dir data/CellPose \
  --image-ext .tif \
  --scale 4
```

## Import

```python
from src.segmentation_simple import segment_image
```
