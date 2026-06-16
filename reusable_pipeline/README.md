# Cachectic spatial pipeline

This folder contains a small reusable pipeline for processing Visium HD samples.
Nuclei and fibres are handled by the same segmented-object workflow.

Edit `src/cachectic_spatial/config.py` before running the example.

```bash
cd reusable_pipeline
pip install -e .
python3 run_pipeline.py
```

The workflow adds:

- tissue annotations and fluorescence values to Visium HD bin tables;
- labels, polygons, transcript counts, fluorescence values, and morphology for nuclei;
- the same information for fibres.
