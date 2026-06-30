# Mask To Polygons

Convert a labeled TIFF mask into GeoJSON polygons with a `label_id` column.

## CLI

```bash
python -m src.polygons \
  --image-mask data/ARVIS/Project_BLOCCO_1_ch00.tif \
  --output data/ARVIS/Project_BLOCCO_1_ch00_polygons.geojson
```

## Import

```python
from src.polygons import to_polygons, save_gdf
```
