# Tools

Small command-line tools and importable Python modules used in this project.

Each tool directory contains:

- `__main__.py` for `python -m src.<tool_or_package>`
- `__init__.py` for Python imports
- `README.md` with minimal usage
- `requirements.txt` with only the packages needed by that tool

Current tools:

- `background`: create a tissue/background mask from an RGB TIFF.
- `segmentation`: multiscale Cellpose segmentation with tile rescue.
- `segmentation_simple`: simple Cellpose segmentation wrapper.
- `polygons`: convert labeled TIFF masks to GeoJSON polygons.
- `utils`: split/merge TIFF and merge GeoJSON using sample bounding boxes.
- `spatial_data.build`: build block/sample SpatialData objects.
- `spatial_data.merge`: merge saved SpatialData zarr stores.
- `spatial_data.qc`: plot SpatialData or raw raster QC.

Typical block workflow:

1. Merge sample-level ARVIS TIFFs with `python -m src.utils tiff merge`.
2. Merge sample-level ROI GeoJSON files with `python -m src.utils geojson merge`.
3. Run or collect Cellpose masks.
4. Convert ARVIS and Cellpose masks to polygons with `python -m src.polygons`.
5. Create background polygons with `python -m src.background`.
6. Build block and sample SpatialData zarr stores with `python -m src.spatial_data.build`.
7. Merge multiple sample/block zarr stores with `python -m src.spatial_data.merge` when needed.

Sample files are GeoJSON FeatureCollections compatible with Napari polygon exports. Each sample is a polygon feature and must include the block name in `properties.block`:

```json
"type": "Feature",
"geometry": {
  "type": "Polygon",
  "coordinates": [[[0, 0], [16166, 0], [16166, 8000], [0, 8000]]]
},
"properties": {
  "object_type": "annotation",
  "isLocked": false,
  "block": "blocco1",
  "sample": "c26STAT3",
  "name": "c26STAT3",
  "sample_key": "blocco1_c26STAT3"
}
```

Current scripts use the external bounding box of each polygon for file crops and SpatialData queries, and preserve the polygon-derived rotated rectangle in sample `attrs`.
