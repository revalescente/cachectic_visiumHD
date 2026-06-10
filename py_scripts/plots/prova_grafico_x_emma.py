import geopandas as gpd
import pandas as pd
import numpy as np
import spatialdata as sd
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa

# Grafico con i nuovi dati

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/arivis_plus_bins/blocco4_c26")

# mi servono:
#             nuclei Trim63
#             bin IIb azzurri
#             bin myo grigi

# prova 1, creo 2 colonne T/F 
# prova 2, creo 2 colonne, una sulla tabella dei nuclei con trim63 1/0, e una su tab bin 1,2,0

shape_nuclei = "nuclei_arivis_poly"
shape_bin = "blocco4_square_016um"

table_bin = "square_016um"
table_nuclei = "arivis_nuclei_table"
# create the first column
# sdata["arivis_nuclei_table"].obs["trim_check"] = sdata["arivis_nuclei_table"].obs["nuclei_types"] == "Myonuclei_Trim63"
# sdata["square_016um"].obs["IIb_check"] = sdata["square_016um"].obs["bin_types"] == "Myonuclei_IIb"

# Creates a column of 1s and 0s
sdata[table_nuclei].obs["trim_check"] = (sdata[table_nuclei].obs["nuclei_types"] == "Myonuclei_Trim63").astype(int)
sdata[table_bin].obs["IIb_check"] = (sdata[table_bin].obs["bin_types"] == "Myonuclei_IIb").astype(int)
sdata[table_bin].obs["myo_check"] = (sdata[table_bin].obs["bin_types"] != "Myonuclei_IIb").astype(int)

# to have invisible shapes where the values == 0 use this cmap
#user_def_cmap = "cividis"
nuclei_color = "plasma"
cmap_nuclei = set_zero_in_cmap_to_transparent(cmap=nuclei_color)
IIb_color = "cool"
cmap_IIb = set_zero_in_cmap_to_transparent(cmap=IIb_color)
myo_color = "grey"
cmap_myo = set_zero_in_cmap_to_transparent(cmap="grey")

# 2. Filter for IIb_check == 1
table_filtered = sdata[table_bin][sdata[table_bin].obs["IIb_check"] == 1].copy()

sdata["table_IIb"] = table_filtered
prova = sd.match_element_to_table(sdata, shape_bin, "table_IIb")
sdata['IIb_to_plot'] = prova[0]['blocco4_square_016um']

sdata['table_IIb'].uns = {'spatialdata_attrs': {'instance_key': 'location_id', 'region': ['IIb_to_plot'], 'region_key': 'region'}}
sdata['table_IIb'].obs['region'] = "IIb_to_plot"
sdata['table_IIb'].obs['region'] = sdata['table_IIb'].obs['region'].astype('category')


# final plot
fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_facecolor('black')        # Set the plot area to black

sdata.pl.render_images(
  "blocco4_full_image"
).pl.render_shapes(
    'IIb_to_plot', 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="lightblue",
    #cmap=cmap_IIb,
    table_name = "table_IIb"
).pl.render_shapes(
    shape_nuclei, 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="trim_check",
    cmap = cmap_nuclei
).pl.show(
    ax=ax, 
    title="Seconda prova grafico emma", 
    coordinate_systems="blocco4" # make sure this matches your coordinate system
)

# 6. Save and close
save_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_forsefinale_img.png"
fig.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)




# add grey mionuclei
mask = (
    sdata[table_bin].obs["bin_types"].str.startswith("Myonuclei_") & 
    (sdata[table_bin].obs["bin_types"] != "Myonuclei_IIb")
)
# Convert True/False to 1/0 and assign to a new column
sdata[table_bin].obs["Myo_other_check"] = mask.astype(int)

# Verify the result
print(sdata[table_bin].obs.groupby("bin_types")["Myo_other_check"].value_counts())

table_filtered = sdata[table_bin][sdata[table_bin].obs["Myo_other_check"] == 1].copy()

sdata["table_myo"] = table_filtered
prova = sd.match_element_to_table(sdata, shape_bin, "table_myo")
sdata['Myo_to_plot'] = prova[0]['blocco4_square_016um']

sdata['table_myo'].uns = {'spatialdata_attrs': {'instance_key': 'location_id', 'region': ['Myo_to_plot'], 'region_key': 'region'}}
sdata['table_myo'].obs['region'] = "Myo_to_plot"
sdata['table_myo'].obs['region'] = sdata['table_myo'].obs['region'].astype('category')


# final plot
fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_images(
  "blocco4_full_image"
).pl.render_shapes(
    'IIb_to_plot', 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="lightblue",
    #cmap=cmap_IIb,
    table_name = "table_IIb"
).pl.render_shapes(
    shape_nuclei, 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="trim_check",
    cmap = cmap_nuclei
).pl.show(
    ax=ax, 
    title="Grafico per emma", 
    coordinate_systems="blocco4" # make sure this matches your coordinate system
)
# 6. Save and close
save_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_emma_conimg.png"
fig.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)

# final plot
fig, ax = plt.subplots(figsize=(20, 20))
#ax.set_facecolor('black')        # Set the plot area to black
sdata.pl.render_shapes(
    'IIb_to_plot', 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="lightblue",
    #cmap=cmap_IIb,
    table_name = "table_IIb"
).pl.render_shapes(
    'Myo_to_plot', 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="grey",
    #cmap=cmap_IIb,
    table_name = "table_myo"
).pl.render_shapes(
    shape_nuclei, 
    outline=False, 
    outline_alpha=0, 
    fill_alpha=1,
    color="trim_check",
    cmap = cmap_nuclei
).pl.show(
    ax=ax, 
    title="Grafico per emma", 
    coordinate_systems="blocco4" # make sure this matches your coordinate system
)
# 6. Save and close
save_path = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/grafico_emma_noimg.png"
fig.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
