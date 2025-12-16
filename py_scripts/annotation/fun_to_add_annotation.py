import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import pandas as pd
import argparse
from pathlib import Path
import pandas as pd
import spatialdata as sd


def add_annotations_to_nuclei_counts(sdata_path: str, annotations_path: str, save: bool = True):
    """
    Add annotation results to the nuclei_counts table in a SpatialData object.
    
    Parameters
    ----------
    sdata_path : str
        Path to the SpatialData . zarr store. 
    annotations_path : str
        Path to the annotation results .parquet file. 
    save : bool
        Whether to save the updated SpatialData object back to disk.
    
    Returns
    -------
    sdata : SpatialData
        The updated SpatialData object. 
    """
    # Load SpatialData
    print(f"Loading SpatialData from {sdata_path}...")
    sdata = sd.read_zarr(sdata_path)
    
    # Load annotation results
    print(f"Loading annotations from {annotations_path}...")
    annotations = pd.read_parquet(annotations_path)
    
    # Set cell_id as index for annotations
    if 'cell_id' in annotations.columns:
        annotations = annotations.set_index('cell_id')
    elif 'row_id' in annotations.columns:
        # In case you saved row names as 'row_id' during R export
        annotations = annotations.set_index('row_id')
    else:
        raise KeyError("Annotation file must contain 'cell_id' or 'row_id' column.")
    
    print(f"Annotations shape: {annotations.shape}")
    
    # Get the nuclei_counts table
    nuclei_counts = sdata.tables['nuclei_counts']
    print(f"nuclei_counts shape: {nuclei_counts.shape}")
    
    # Ensure nuclei_counts has cell_id as index in . obs
    # If cell_id is a column, set it as index
    if 'cell_id' in nuclei_counts.obs.columns:
        nuclei_counts.obs = nuclei_counts.obs.set_index('cell_id')
    
    print(f"nuclei_counts obs index: {nuclei_counts.obs.index[:5]. tolist()}...")
    print(f"annotations index: {annotations.index[:5].tolist()}...")
    
    # Find common indices
    common_ids = nuclei_counts.obs.index.intersection(annotations.index)
    print(f"Common cell_ids: {len(common_ids)} / {len(nuclei_counts.obs.index)}")
    
    if len(common_ids) == 0:
        raise ValueError("No matching cell_ids found between nuclei_counts and annotations!")
    
    # Get annotation columns (exclude any that already exist in . obs to avoid duplicates)
    existing_cols = set(nuclei_counts.obs.columns)
    new_cols = [col for col in annotations.columns if col not in existing_cols]
    print(f"Adding {len(new_cols)} new annotation columns: {new_cols[:10]}{'...' if len(new_cols) > 10 else ''}")
    
    # Merge annotations into . obs
    # Use . loc to align by index and handle missing values
    # for col in new_cols:
    #     nuclei_counts.obs[col] = annotations. loc[nuclei_counts.obs.index. intersection(annotations.index), col].reindex(nuclei_counts.obs.index)
    nuclei_counts.obs = nuclei_counts.obs.join(annotations[new_cols], how='left')
    
    print(f"Updated nuclei_counts. obs shape: {nuclei_counts.obs.shape}")
    
    # Update the table in sdata
    sdata.tables['nuclei_counts'] = nuclei_counts
    
    # Save the updated SpatialData object
    if save:
        print(f"Saving updated SpatialData to {sdata_path}...")
        # Overwrite the existing table
        sdata.delete_element_from_disk('nuclei_counts')
        sdata.write_element('nuclei_counts', overwrite=True)
        print("Done!")
    
    return sdata


def batch_add_annotations(sdata_dir: str, annotations_dir: str):
    """
    Batch process multiple SpatialData stores and their corresponding annotation files.
    
    Assumes matching names: sample1. zarr <-> sample1.parquet
    """
    sdata_path = Path(sdata_dir)
    annot_path = Path(annotations_dir)
    
    for zarr_store in sdata_path.glob("*.zarr"):
        sample_name = zarr_store.stem  # e.g., "blocco1_sham"
        parquet_file = annot_path / f"{sample_name}.parquet"
        
        if parquet_file.exists():
            print(f"\n{'='*60}")
            print(f"Processing {sample_name}")
            print(f"{'='*60}")
            try:
                add_annotations_to_nuclei_counts(str(zarr_store), str(parquet_file))
            except Exception as e:
                print(f"Error processing {sample_name}: {e}")
        else:
            print(f"Warning: No annotation file found for {sample_name} (expected {parquet_file})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add annotation results to SpatialData nuclei_counts table.")
    parser.add_argument("--sdata", type=str, help="Path to a single SpatialData . zarr store")
    parser.add_argument("--annotations", type=str, help="Path to a single annotations .parquet file")
    parser.add_argument("--sdata_dir", type=str, help="Path to directory containing multiple . zarr stores (for batch processing)")
    parser.add_argument("--annotations_dir", type=str, help="Path to directory containing multiple .parquet files (for batch processing)")
    
    args = parser.parse_args()
    
    if args.sdata and args.annotations:
        # Single file mode
        add_annotations_to_nuclei_counts(args.sdata, args.annotations)
    elif args.sdata_dir and args.annotations_dir:
        # Batch mode
        batch_add_annotations(args.sdata_dir, args.annotations_dir)
    else:
        print("Usage:")
        print("  Single file:  python add_annotations_to_sdata.py --sdata /path/to/sample.zarr --annotations /path/to/sample.parquet")
        print("  Batch mode:   python py_scripts/annotation/fun_to_add_annotation.py --sdata_dir /mnt/europa/valerio/data/zarr_store/samples --annotations_dir /mnt/europa/valerio/data/data_tables/annotation_results_filtered")
