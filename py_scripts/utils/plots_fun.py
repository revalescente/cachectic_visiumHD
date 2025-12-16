import os
import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot
import pandas as pd

def sample_cells_by_type(
    sdata,
    table_key,
    which_cells,
    label_col='labels_singler',
    n_samples=5,
    random_state=42
):
    """
    Sample cells from specific cell types in a SpatialData object. 

    Parameters
    ----------
    sdata : SpatialData
        SpatialData object containing the table.
    table_key : str
        Key for the table in sdata.tables.
    which_cells : list
        List of cell types to sample from.
    label_col :  str, optional
        Column name containing cell type labels (default: 'labels_singler').
    n_samples : int, optional
        Maximum number of cells to sample per cell type (default: 5).
    random_state : int, optional
        Random seed for reproducibility (default: 42).

    Returns
    -------
    pd.DataFrame
        DataFrame containing sampled cells from all specified cell types.
    """
    # Access the obs DataFrame
    obs_df = sdata.tables[table_key].obs

    # Initialize list to store sampled subsets
    collected_samples = []

    # Loop through cell types and sample
    for cell_type in which_cells:
        # Filter for the specific cell type
        subset = obs_df[obs_df[label_col] == cell_type]
        n_available = len(subset)

        if n_available > 0:
            # Sample up to n_samples
            n_sample = min(n_samples, n_available)
            sampled_subset = subset.sample(n=n_sample, random_state=random_state)
            collected_samples.append(sampled_subset)
            print(f"Sampled {n_sample}/{n_available} cells for '{cell_type}'")
        else:
            print(f"Warning: No cells found for '{cell_type}'")

    # Concatenate into a single DataFrame
    if collected_samples:
        summary_df = pd.concat(collected_samples)
        print(f"\nTotal sampled:  {len(summary_df)} cells")
        return summary_df
    else: 
        print("No cells sampled.")
        return pd.DataFrame()

def plot_single_nuclei(
    df,
    sdata,
    output_dir,
    image_key,
    shapes_key,
    table_key,
    coord_sys,
    buffer=50,
    figsize=(10, 10),
    cell_id_col='cell_id',
    cell_type_col='labels_singler',
    color_key = "labels_singler",
    highlight_color='yellow',
    highlight_width=2.5,
    outline_alpha=0.5,
    outline_width=1.5,
    fill_alpha=1
):
    """
    Plot individual cells from a dataframe with their surrounding context. 

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing cell_id and cell_type columns for cells to plot.
    sdata : SpatialData
        SpatialData object containing shapes, images, and tables.
    output_dir : str
        Directory to save the output PNG files.
    shapes_key :  str
        Key for the nuclei/cell shapes in sdata.shapes.
    image_key : str
        Key for the image to render in sdata.images.
    coord_sys : str
        Coordinate system to use for plotting. 
    buffer : float, optional
        Buffer around cell bounding box in pixels (default: 50).
    figsize : tuple, optional
        Figure size for plots (default: (10, 10)).
    cell_id_col : str, optional
        Column name for cell IDs in df (default: 'cell_id').
    cell_type_col : str, optional
        Column name for cell types in df (default: 'labels_singler').
    highlight_color :  str, optional
        Color for highlighting the target cell (default: 'yellow').
    highlight_width : float, optional
        Outline width for highlighted cell (default: 2).
    outline_alpha : float, optional
        Alpha for background cell outlines (default: 0.5).
    outline_width : float, optional
        Outline width for background cells (default: 1.0).

    Returns
    -------
    dict
        Summary with 'success' and 'failed' lists of cell_ids.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    temp_shape_key = "_temp_highlight_shape"
    results = {'success': [], 'failed': []}

    # Loop through the cells
    for index, row in df.iterrows():
        cell_id = row[cell_id_col]
        cell_type = row[cell_type_col]
        save_path = os.path.join(output_dir, f"{cell_type}_{cell_id}.png")
        print(f"Processing {cell_id} ({cell_type})...")

        try:
            # 1. Check if cell_id exists in shapes
            if cell_id not in sdata.shapes[shapes_key].index:
                raise KeyError(f"Cell ID {cell_id} not found in shapes")

            # 2. Isolate the specific shape for this cell
            sdata.shapes[temp_shape_key] = sdata.shapes[shapes_key].loc[[cell_id]].copy()

            # 3. Calculate the extent (bounding box) of this single cell
            extent = sd.get_extent(sdata.shapes[temp_shape_key], coordinate_system=coord_sys)
            min_x, max_x = extent['x']
            min_y, max_y = extent['y']

            # 4. Define the plot window (cell bbox + buffer)
            bbox_min = [min_x - buffer, min_y - buffer]
            bbox_max = [max_x + buffer, max_y + buffer]

            # 5. Setup plot
            plt.figure(figsize=figsize)
            ax = plt.gca()

            # 6. Query, render, and save
            sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=bbox_min,
                max_coordinate=bbox_max,
                target_coordinate_system=coord_sys,
            ).pl.render_images(
                image_key
            ).pl.render_shapes(
                shapes_key,
                outline=True,
                outline_alpha=outline_alpha,
                outline_width=outline_width,
                fill_alpha=fill_alpha,
                color = color_key,
                table_name = table_key
            ).pl.render_shapes(
                temp_shape_key,
                fill_alpha=0,
                outline=True,
                outline_color=highlight_color,
                outline_width=highlight_width,
                outline_alpha=1
            ).pl.show(
                ax=ax,
                title=f"{cell_type}\n{cell_id}",
                coordinate_systems=coord_sys,
                save=save_path
            )

            plt.close()
            results['success'].append(cell_id)

        except KeyError as e: 
            print(f"  Error: {e}")
            results['failed'].append(cell_id)
        except Exception as e:
            print(f"  Error plotting {cell_id}: {e}")
            results['failed'].append(cell_id)
        finally:
            # Cleanup the temporary shape from sdata
            if temp_shape_key in sdata.shapes:
                del sdata.shapes[temp_shape_key]
            plt.close('all')

    # Summary
    print(f"\nCompleted:  {len(results['success'])} success, {len(results['failed'])} failed")
    return results
