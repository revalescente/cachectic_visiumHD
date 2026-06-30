# SpatialData Build

Build a block-level SpatialData object from SpaceRanger HD outputs, merged ARVIS/Cellpose polygons, ROI annotations, and background polygons. The tool also saves sample-level crops using `samples_he.json`.

`--wga-rgb-image` is intended for `data/WGA_GFP_DAPI/RGB/Project_BLOCCO_<n>.tif`. Until those images are aligned to H&E, pass the H&E image itself as a placeholder for pipeline tests.

Sample files can be Napari-style GeoJSON FeatureCollections. The build uses each sample polygon's external bounding box for the crop and stores the polygon-derived rotated rectangle in each sample SpatialData `attrs`.

## CLI

```bash
python -m src.spatial_data.build blocco1 \
  --samples-json data/samples_he.json \
  --spaceranger data/SpaceRanger/blocco1/outs \
  --fullres-image data/H\&E/Project_BLOCCO_1_ch00.tif \
  --wga-rgb-image data/H\&E/Project_BLOCCO_1_ch00.tif \
  --nuclei-polygons data/ARVIS/Project_BLOCCO_1_ch00_polygons.geojson \
  --fibres-polygons data/CellPose/Project_BLOCCO_1_ch00_polygons.geojson \
  --background-polygons data/H\&E/Project_BLOCCO_1_background_mask.geojson \
  --roi-metadata data/ARVIS/Project_BLOCCO_1_ch00.geojson \
  --output-dir data/spatialdata
```

## Import

```python
from src.spatial_data.build import build_block, build_block_and_samples, crop_samples
```
