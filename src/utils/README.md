# Split/Merge

Split and merge block/sample files using sample bounding boxes from a JSON file.

Sample files can be Napari-style GeoJSON FeatureCollections. Each feature must contain `properties.block`; `properties.sample`, `properties.name`, and `properties.sample_key` are used when present. The tools use each polygon's external bounding box. The legacy `min_coordinate`/`max_coordinate` dictionary format remains supported.

Supported data types:

- TIFF: split block images into sample images; merge sample images back into block coordinates.
- GeoJSON: merge sample-level annotations into block coordinates.

## CLI

```bash
python -m src.utils tiff split blocco1 \
  --json data/samples_he.json \
  --image data/H\&E/Project_BLOCCO_1_ch00.tif \
  --output-dir data/H\&E/Project_BLOCCO_1_ch00
```

```bash
python -m src.utils geojson merge blocco1 \
  --json data/samples_he.json \
  --geojson-dir data/ARVIS/intissue_GFP_blocco1 \
  --output data/ARVIS/Project_BLOCCO_1_ch00.geojson
```

## Import

```python
from src.utils import split_tiff, merge_tiff, merge_geojson
```
