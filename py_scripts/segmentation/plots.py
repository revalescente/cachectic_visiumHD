

def plot_nuclei(
    sdata: sd.SpatialData,
    extreme_nuclei_df: pd.DataFrame,
    shapes_key: str,
    save_dir: str = 'output_python/nuclei_Strange/',
    padding: int = 50,
    figsize: tuple = (15, 15)
):
    """
    Plots individual nuclei from a summary DataFrame and saves them as descriptive images

    Args:
        sdata: The SpatialData object.
        extreme_nuclei_df: DataFrame with 'cell_id', 'feature', 'extreme_type'.
        shapes_key: The key for the nuclei shapes (e.g., 'blocco4_nuclei_boundaries').
        save_dir: The base directory where plot subdirectories will be created.
        padding: Pixels to add around the nucleus bounding box.
        figsize: The size of the output figure.
    """
    # 1. Get the full GeoDataFrame of all nuclei
    all_nuclei_gdf = sdata[shapes_key]
    
    # 2. Infer blocco_key and image_key from shapes_key
    match = re.search(r'(blocco\d+)', shapes_key)
    if not match:
        raise ValueError(f"Could not infer 'blocco_key' from shapes_key '{shapes_key}'.")
    blocco_key = match.group(1)
    image_key = f"{blocco_key}_full_image"
    print(f"Inferred blocco_key: '{blocco_key}' and image_key: '{image_key}'")

    # 3. Loop through each row of the DataFrame
    for index, row in extreme_nuclei_df.iterrows():
        cell_id = row['cell_id']
        feature = row['feature']
        extreme_type = row['extreme_type'].replace(' ', '')

        # Create subdirectory for the feature
        feature_save_dir = Path(save_dir) / feature
        feature_save_dir.mkdir(parents=True, exist_ok=True)
        
        # Construct descriptive filename
        save_filename = f"{feature}_{extreme_type}_{cell_id}.png"
        save_path = feature_save_dir / save_filename

        print(f"Processing nucleus: {cell_id} for feature: {feature}...")

        try:
            # 4. Select the single nucleus
            single_nucleus_gdf = all_nuclei_gdf.loc[[cell_id]]

            # --- YOUR PLOTTING DYNAMIC ---
            # Add the specific nucleus as a temporary shape to sdata
            sdata['nucleo_shape'] = single_nucleus_gdf

            # 5. Calculate its bounding box and add padding
            extent = sd.get_extent(sdata['nucleo_shape'], coordinate_system=blocco_key, exact=True)
            min_coord = [extent['x'][0] - padding, extent['y'][0] - padding]
            max_coord = [extent['x'][1] + padding, extent['y'][1] + padding]

            # 6. Create figure and axis
            fig, ax = plt.subplots(figsize=figsize)
            
            # 7. Query, render, and show using your exact command structure
            sdata.query.bounding_box(
                axes=["x", "y"],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=blocco_key,
            ).pl.render_images(image_key).pl.render_shapes(
                shapes_key, outline=True, outline_width=1.5, fill_alpha=0
            ).pl.render_shapes(
              'nucleo_shape', outline=True, color="red", fill_alpha=0.5
            ).pl.show(ax=ax, coordinate_systems=blocco_key, save=save_path)
            
            plt.close(fig) # Close the figure to free up memory

        except KeyError:
            print(f"Warning: cell_id '{cell_id}' not found in '{shapes_key}'. Skipping.")
        except Exception as e:
            print(f"Error plotting '{cell_id}': {e}")
        finally:
            # --- CRITICAL CLEANUP STEP ---
            # Remove the temporary shape to ensure sdata is clean for the next loop
            if 'nucleo_shape' in sdata:
                del sdata['nucleo_shape']
            
    print("\nDone.")
