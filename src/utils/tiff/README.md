# TIFF Split/Merge

Split a block-level TIFF into sample TIFFs, or merge sample TIFFs back into block coordinates using a samples JSON file.

## CLI

```bash
python -m src.utils tiff split blocco1 \
  --json data/samples_he.json \
  --image data/H\&E/Project_BLOCCO_1_ch00.tif \
  --output-dir data/H\&E/Project_BLOCCO_1_ch00
```

```bash
python -m src.utils tiff merge blocco1 \
  --json data/samples_he.json \
  --images-dir data/H\&E/Project_BLOCCO_1_ch00 \
  --output data/H\&E/Project_BLOCCO_1_ch00_rec.tif
```

## Import

```python
from src.utils.tiff import split, merge
```
