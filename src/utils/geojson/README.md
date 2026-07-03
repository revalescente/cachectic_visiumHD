# GeoJSON Split/Merge

Merge sample-level GeoJSON annotations into block-level coordinates using sample bounding boxes from a JSON file.

## CLI

```bash
python -m src.utils geojson merge blocco1 \
  --json data/samples_he.json \
  --geojson-dir data/ARVIS/intissue_GFP_blocco1 \
  --output data/ARVIS/Project_BLOCCO_1_ch00.geojson
```

## Import

```python
from src.utils.geojson import merge
```
