# SpatialData QC

Plot quick QC figures from a saved SpatialData zarr store. A raw TIFF mode is also available for pre-build checks.

## CLI

```bash
python -m src.spatial_data.qc \
  --zarr data/spatialdata/blocco1_sham.zarr \
  --coordinate-system blocco1_sham
```

Raw pre-build QC:

```bash
python -m src.spatial_data.qc blocco1 \
  --samples-json data/samples_he.json \
  --he data/H\&E/Project_BLOCCO_1_ch00.tif \
  --nuclei-labels data/ARVIS/Project_BLOCCO_1_ch00.tif \
  --fibres-labels data/CellPose/Project_BLOCCO_1_ch00.tif
```

## Import

```python
from src.spatial_data.qc import plot_zarr, summarize_zarr
```
