# SpatialData Merge

Merge multiple saved SpatialData zarr stores into a single object with prefixed element keys.

## CLI

```bash
python -m src.spatial_data.merge \
  --output data/spatialdata/samples_merged.zarr \
  blocco1_sham=data/spatialdata/blocco1_sham.zarr \
  blocco2_c26=data/spatialdata/blocco2_c26.zarr
```

## Import

```python
from src.spatial_data.merge import merge_sdata_objects
```
