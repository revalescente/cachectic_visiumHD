# make_plots.py
from shapes_plots_utils import custom_spatial_plot

# BLOCCO 4 ----------------------------------------------------------------
#### PLOT 1 - 1 cell type bin - 1 nucleus cell type - with image

my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "lightblue" # Use "color" for solid matplotlib colors
    },

    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "plasma" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco4_c26",
    "shape_bin": "blocco4_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco4"
}

# Plot 
custom_spatial_plot(
    **common_args, # Unpacks all the arguments above
    image_name="blocco4_full_image",
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_con_img.png"
)


#### PLOT 2 - 3 cell types bin - 1 nucleus cell type - no image


my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "blue" # Use "color" for solid matplotlib colors
    },
    {
        "name": "Bin_Other_Myo", 
        "types_startswith": "Myonuclei_", 
        "exclude_types": ["Myonuclei_IIb", "Myonuclei_Trim63"], 
        "color": "lightgrey"
    },
    {
        "name": "Bin_trim", # free choice
        "types": ["Myonuclei_Trim63"], # from the bin_types column
        "color": "orange" # Use "color" for solid matplotlib colors
    }
    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "OrRd" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco4_c26",
    "shape_bin": "blocco4_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco4"
}

# Plot 2
custom_spatial_plot(
    **common_args,
    image_name=None, # Set to None to disable image
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_senza_img_3bins_1nuc.png"
)

# BLOCCO 6 ----------------------------------------------------------------


my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "blue" # Use "color" for solid matplotlib colors
    },
    {
        "name": "Bin_Other_Myo", 
        "types_startswith": "Myonuclei_", 
        "exclude_types": ["Myonuclei_IIb", "Myonuclei_Trim63"], 
        "color": "lightgrey"
    },
    {
        "name": "Bin_trim", # free choice
        "types": ["Myonuclei_Trim63"], # from the bin_types column
        "color": "orange" # Use "color" for solid matplotlib colors
    }
    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "OrRd" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco6_c26",
    "shape_bin": "blocco6_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco6"
}

# Plot
custom_spatial_plot(
    **common_args,
    image_name=None, # Set to None to disable image
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/blocco6_senza_img_3bins_1nuc.png"
)

#### PLOT 2 - 2 cell types bin - 1 nucleus cell type - no image

my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "lightblue" # Use "color" for solid matplotlib colors
    },

    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "plasma" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco6_c26",
    "shape_bin": "blocco6_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco6"
}

# Plot 
custom_spatial_plot(
    **common_args, # Unpacks all the arguments above
    image_name="blocco6_full_image",
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/blocco6_con_img.png"
)


# BLOCCO 2 ----------------------------------------------------------------

my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "blue" # Use "color" for solid matplotlib colors
    },
    {
        "name": "Bin_Other_Myo", 
        "types_startswith": "Myonuclei_", 
        "exclude_types": ["Myonuclei_IIb", "Myonuclei_Trim63"], 
        "color": "lightgrey"
    },
    {
        "name": "Bin_trim", # free choice
        "types": ["Myonuclei_Trim63"], # from the bin_types column
        "color": "orange" # Use "color" for solid matplotlib colors
    }
    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "OrRd" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco2_c26",
    "shape_bin": "blocco2_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco2"
}

# Plot
custom_spatial_plot(
    **common_args,
    image_name=None, # Set to None to disable image
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/blocco2_senza_img_3bins_1nuc.png"
)

#### PLOT 2 - 2 cell types bin - 1 nucleus cell type - no image

my_bins = [
    {
        "name": "Bin_IIb", # free choice
        "types": ["Myonuclei_IIb"], # from the bin_types column
        "color": "lightblue" # Use "color" for solid matplotlib colors
    },

    # Want to add Endothelial? Just add another block:
    # { "name": "Bin_Endo", "types": ["Endothelial"], "color": "yellow" }
]

my_nuclei = [
    {
        "name": "Nuclei_Trim63", 
        "types": ["Myonuclei_Trim63"], 
        "cmap": "plasma" # Use "cmap" (colormap) so 0 values become transparent
    }
    # Want to add another nuclei type?
    # { "name": "Nuclei_FAP", "types": ["FAPs"], "cmap": "cool" }
]

common_args = {
    "sdata_path": "/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco2_c26",
    "shape_bin": "blocco2_square_016um",
    "table_bin": "square_016um",
    "bin_obs_col": "bin_types",
    "bin_configs": my_bins,
    "shape_nuclei": "nuclei_arivis_poly",
    "table_nuclei": "arivis_nuclei_table",
    "nuclei_obs_col": "nuclei_types",
    "nuclei_configs": my_nuclei,
    "coordinate_system": "blocco2"
}

# Plot 
custom_spatial_plot(
    **common_args, # Unpacks all the arguments above
    image_name="blocco2_full_image",
    save_path="/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/blocco2_con_img.png"
)
