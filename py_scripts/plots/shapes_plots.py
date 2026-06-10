import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import spatialdata_plot
import matplotlib.pyplot as plt


def create_filtered_element(sdata, original_shape, original_table, mask, element_name):
    """
    Filters a table based on a mask, subsets the corresponding shapes, 
    and wires the metadata correctly in the SpatialData object.
    """
    table_filtered = sdata[original_table][mask].copy()
    new_table_name = f"table_{element_name}"
    sdata[new_table_name] = table_filtered
    
    # Dynamically get the original instance key instead of hardcoding 'location_id'
    original_attrs = sdata[original_table].uns.get('spatialdata_attrs', {})
    instance_key = original_attrs.get('instance_key', 'location_id') # fallback
    
    # Match the filtered table back to the original shapes
    matched = sd.match_element_to_table(sdata, original_shape, new_table_name)
    sdata[element_name] = matched[0][original_shape]

    # Fix the spatialdata_attrs metadata with the correct instance_key
    sdata[new_table_name].uns['spatialdata_attrs'] = {
        'instance_key': instance_key, 
        'region': [element_name], 
        'region_key': 'region'
    }
    sdata[new_table_name].obs['region'] = element_name
    sdata[new_table_name].obs['region'] = sdata[new_table_name].obs['region'].astype('category')


def custom_spatial_plot(
    sdata_path,
    shape_bin, table_bin, bin_obs_col, bin_configs,
    shape_nuclei, table_nuclei, nuclei_obs_col, nuclei_configs,
    coordinate_system="blocco4",
    image_name=None,      # Set to None if you don't want an image
    save_path=None
):
    print(f"Loading data from {sdata_path}...")
    sdata = sd.read_zarr(sdata_path)
    
    fig, ax = plt.subplots(figsize=(20, 20))
    
    # ==========================================
    # PHASE 1: PREPARE ALL DATA
    # ==========================================
    # 1A. Prepare Bins (using the filtering method)
    for config in bin_configs:
        name = config["name"]
        print(f"Preparing bin group: {name}")
        
        if "types" in config:
            mask = sdata[table_bin].obs[bin_obs_col].isin(config["types"])
        elif "types_startswith" in config:
            mask = sdata[table_bin].obs[bin_obs_col].str.startswith(config["types_startswith"])
            if "exclude_types" in config:
                mask = mask & (~sdata[table_bin].obs[bin_obs_col].isin(config["exclude_types"]))
        else:
            continue
            
        create_filtered_element(sdata, shape_bin, table_bin, mask, name)
        if config.get("cmap"):
            sdata[f"table_{name}"].obs["plot_val"] = 1

    # 1B. Prepare Nuclei (using YOUR 1/0 binary column method)
    for config in nuclei_configs:
        name = config["name"]
        print(f"Preparing nuclei group: {name} (1/0 mask)")
        
        # Create a column of 1s and 0s directly on the original nuclei table
        mask = sdata[table_nuclei].obs[nuclei_obs_col].isin(config["types"])
        sdata[table_nuclei].obs[f"{name}_check"] = mask.astype(int)

    # ==========================================
    # PHASE 2: BUILD THE PLOTTING TREE
    # ==========================================
    if image_name is not None:
        plot_obj = sdata.pl.render_images(image_name)
    else:
        plot_obj = sdata

    # 2A. Render Bins
    for config in bin_configs:
        name = config["name"]
        kwargs = {"outline": False, "outline_alpha": 0, "fill_alpha": 1, "table_name": f"table_{name}"}
        
        if config.get("cmap"):
            kwargs["color"] = "plot_val"
            kwargs["cmap"] = set_zero_in_cmap_to_transparent(cmap=config["cmap"])
        else:
            kwargs["color"] = config["color"]
            
        plot_obj = plot_obj.pl.render_shapes(name, **kwargs)

    # 2B. Render Nuclei
    for config in nuclei_configs:
        name = config["name"]
        # Use the original shape_nuclei, not a filtered one!
        kwargs = {"outline": False, "outline_alpha": 0, "fill_alpha": 1}
        
        # Color by the new 1/0 column we just made
        kwargs["color"] = f"{name}_check"
        kwargs["cmap"] = set_zero_in_cmap_to_transparent(cmap=config.get("cmap", "plasma"))
            
        plot_obj = plot_obj.pl.render_shapes(shape_nuclei, **kwargs)

    # ==========================================
    # PHASE 3: SHOW AND SAVE
    # ==========================================
    print("Generating plot...")
    plot_obj.pl.show(ax=ax, title="Custom Spatial Plot", coordinate_systems=coordinate_system)
    
    if save_path:
        print(f"Saving to {save_path}")
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.close(fig)
    print("Done!")

# ==========================================
# USER CONFIGURATION & EXECUTION
# ==========================================
if __name__ == "__main__":
    
    # 1. Define Bins to plot
    # You can mix specific exact 'types' or broad 'types_startswith'
    BIN_CONFIGS = [
        {
            "name": "Bin_IIb", 
            "types": ["Myonuclei_IIb"], 
            "color": "lightblue" # Flat color
        },
        {
            "name": "Bin_Other_Myo", 
            "types_startswith": "Myonuclei_", 
            "exclude_types": ["Myonuclei_IIb"], 
            "color": "grey"      # Flat color
        }
    ]

    # 2. Define Nuclei to plot
    NUCLEI_CONFIGS = [
        {
            "name": "Nuclei_Trim63", 
            "types": ["Myonuclei_Trim63"], 
            "cmap": "plasma"     # Uses colormap (transparent background)
        }
    ]

    # 3. Run script - WITH IMAGE
    custom_spatial_plot(
        sdata_path="/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco4_c26",
        shape_bin="blocco4_square_016um",
        table_bin="square_016um",
        bin_obs_col="bin_types",
        bin_configs=BIN_CONFIGS,
        shape_nuclei="nuclei_arivis_poly",
        table_nuclei="arivis_nuclei_table",
        nuclei_obs_col="nuclei_types",
        nuclei_configs=NUCLEI_CONFIGS,
        image_name="blocco4_full_image", # Include image
        save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_emma_conimg_dynamic.png"
    )

    # 4. Run script - NO IMAGE (Pass None to image_name)
    custom_spatial_plot(
        sdata_path="/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco4_c26",
        shape_bin="blocco4_square_016um",
        table_bin="square_016um",
        bin_obs_col="bin_types",
        bin_configs=BIN_CONFIGS,
        shape_nuclei="nuclei_arivis_poly",
        table_nuclei="arivis_nuclei_table",
        nuclei_obs_col="nuclei_types",
        nuclei_configs=NUCLEI_CONFIGS,
        image_name=None, # No image
        save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_emma_noimg_dynamic.png"
    )
