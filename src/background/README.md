# Background Mask

Create a background mask from an RGB TIFF image and save it as either GeoJSON polygons or a TIFF mask.

## CLI

```bash
python -m src.background \
  --image data/H\&E/Project_BLOCCO_1_ch00.tif \
  --output data/H\&E/Project_BLOCCO_1_background_mask.geojson
```

## Import

```python
from src.background import background_mask, mask_to_polygons
```
