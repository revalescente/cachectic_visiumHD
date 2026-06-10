# spatial_utils.py
import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import spatialdata_plot
import matplotlib.pyplot as plt

def create_filtered_element(sdata, original_shape, original_table, mask, element_name):
    """Filters a table, subsets shapes, and wires metadata."""
    table_filtered = sdata[original_table][mask].copy()
    new_table_name = f"table_{element_name}"
    sdata[new_table_name] = table_filtered
    
    original_attrs = sdata[original_table].uns.get('spatialdata_attrs', {})
    instance_key = original_attrs.get('instance_key', 'location_id')
    
    matched = sd.match_element_to_table(sdata, original_shape, new_table_name)
    sdata[element_name] = matched[0][original_shape]

    sdata[new_table_name].uns['spatialdata_attrs'] = {
        'instance_key': instance_key, 
        'region': [element_name], 
        'region_key': 'region'
    }
    sdata[new_table_name].obs['region'] = element_name
    sdata[new_table_name].obs['region'] = sdata[new_table_name].obs['region'].astype('category')

def custom_spatial_plot(
    sdata_path, shape_bin, table_bin, bin_obs_col, bin_configs,
    shape_nuclei, table_nuclei, nuclei_obs_col, nuclei_configs,
    coordinate_system="blocco4", image_name=None, save_path=None
):
    print(f"Loading data from {sdata_path}...")
    sdata = sd.read_zarr(sdata_path)
    fig, ax = plt.subplots(figsize=(20, 20))
    #ax.set_facecolor('black') # Good for transparent images
    
    # 1A. Prepare Bins
    for config in bin_configs:
        name = config["name"]
        if "types" in config:
            mask = sdata[table_bin].obs[bin_obs_col].isin(config["types"])
        elif "types_startswith" in config:
            mask = sdata[table_bin].obs[bin_obs_col].str.startswith(config["types_startswith"])
            if "exclude_types" in config:
                mask = mask & (~sdata[table_bin].obs[bin_obs_col].isin(config["exclude_types"]))
        else: continue
        create_filtered_element(sdata, shape_bin, table_bin, mask, name)
        if config.get("cmap"): sdata[f"table_{name}"].obs["plot_val"] = 1

    # 1B. Prepare Nuclei (1/0 mask)
    for config in nuclei_configs:
        name = config["name"]
        mask = sdata[table_nuclei].obs[nuclei_obs_col].isin(config["types"])
        sdata[table_nuclei].obs[f"{name}_check"] = mask.astype(int)

    # 2. Build Plotting Tree
    plot_obj = sdata.pl.render_images(image_name) if image_name else sdata

    for config in bin_configs:
        name = config["name"]
        kwargs = {"outline": False, "outline_alpha": 0, "fill_alpha": 1, "table_name": f"table_{name}"}
        if config.get("cmap"):
            kwargs["color"] = "plot_val"
            kwargs["cmap"] = set_zero_in_cmap_to_transparent(cmap=config["cmap"])
        else:
            kwargs["color"] = config["color"]
        plot_obj = plot_obj.pl.render_shapes(name, **kwargs)

    for config in nuclei_configs:
        name = config["name"]
        kwargs = {"outline": False, "outline_alpha": 0, "fill_alpha": 1, "color": f"{name}_check"}
        kwargs["cmap"] = set_zero_in_cmap_to_transparent(cmap=config.get("cmap", "plasma"))
        plot_obj = plot_obj.pl.render_shapes(shape_nuclei, **kwargs)

    # 3. Show and Save
    print("Generating plot...")
    plot_obj.pl.show(ax=ax, title="Custom Spatial Plot", coordinate_systems=coordinate_system)
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved to {save_path}")
