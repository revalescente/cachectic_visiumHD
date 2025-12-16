import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
import geopandas as gpd
import pandas as pd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from py_scripts.utils.utils_fun import read_from_json

# GRAFICAMENTE NON FUNZIONA, MA A LIVELLO OPERATIVO POI SI, QUINDI NON CONTARE SUI GRAFICI DEI POLIGONI PRIMA DEL FILTRAGGIO
# SI FA PRIMA A VERIFICARE DOPO SE CI SONO RIMASTI NUCLEI.. COMUNQUE DA DEFINIRE MEGLIO LE ZONE IN/OUT TISSUE PERCHE CI SONO
# TANTI NUCLEI FUORI.

samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")

blocco = 'blocco1'
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]
# sopa
sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[0]}.zarr")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[1]}.zarr")
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[2]}.zarr")
# spaceranger
sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[0]}")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[1]}")
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[2]}")
# poly
# GFP - inflamed
geojson_path = f"/mnt/europa/valerio/data/json/geojson_dir/GFP_inflamed/{blocco}_GFP_and_inflamed.geojson"
# intissue
geojson_path = f"/mnt/europa/valerio/data/json/geojson_dir/in_tissue/fullres/{blocco}_allpolys_fullres.geojson" # blocco2 solo intissue

gfp_poly = gpd.read_file(geojson_path)
gfp_poly = gfp_poly.set_crs(None, allow_override=True)
# proviamo a separare i polygoni cosi matplot fa meno fatica... ???
#gfp_poly = gfp_poly.explode(index_parts=False).reset_index(drop=True)
gfp_parse = ShapesModel.parse(gfp_poly, transformations={blocco: Identity()})

sdata1.shapes["GFP_poly"] = gfp_parse
sdata2.shapes["GFP_poly"] = gfp_parse
sdata3.shapes["GFP_poly"] = gfp_parse

# Put your objects in a list
sdata_list = [sdata1, sdata2, sdata3]

# plot of the polys
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('GFP_poly', outline=False, outline_alpha=1, outline_width=1, fill_alpha=0
).pl.show(ax = ax, coordinate_systems=blocco, 
    save = f'output_python/roba_da_buttare/{blocco}_gfp_poly.png'
)


# COMBINED PLOT ----

# Initialize the figure and axis
plt.figure(figsize=(20, 20))
ax = plt.gca()

# Iterate through each object and plot it onto the same axis
for sdata_obj in sdata_list:
    # Note: Ensure 'blocco_key' corresponds to a valid coordinate system 
    # inside each specific sdata_obj. If the keys differ per object, 
    # you may need to extract the correct one (e.g., sdata_obj.coordinate_systems[0])
    
    sdata_obj.pl.render_images(
        'full_image', 
        scale='scale3'
    ).pl.render_shapes(
        'GFP_poly', 
        outline=False, 
        outline_alpha=1, 
        outline_width=1, 
        fill_alpha=0
    ).pl.show(
        ax=ax, 
        coordinate_systems=blocco  # Argument for the coordinate system
    )

# Save the combined plot
plt.savefig(f'/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/roba_da_buttare/{blocco}_combined_gfp_poly.png')
plt.show()

# process of cleaning for copilot ----------------------------------------------

# Assuming your objects are in a list

for sdata in sdata_list:
    # 1. Setup: Get table and element references
    table_name = 'nuclei_counts'
    table = sdata[table_name].copy()
    attrs = table.uns['spatialdata_attrs']
    region_name = attrs['region']
    instance_key = attrs['instance_key']
    # Initialize annotation column
    table.obs['in_treatment'] = False
    # Check if polygons exist
    if 'GFP_poly' not in sdata.shapes:
        print(f"Skipping {region_name}: 'GFP_poly' not found.")
        continue
    # 2. Determine Coordinate System (pick the first available)
    nuclei_element = sdata[region_name].copy()
    coord_system = list(get_transformation(nuclei_element, get_all=True).keys())[0]
    # 3. Perform SINGLE Spatial Join
    # Identifies all nuclei overlapping ANY polygon in 'GFP_poly'
    joined_nuclei = sopa.spatial.sjoin(
        sdata, 
        region_name, 
        'GFP_poly', 
        how='inner', 
        predicate='intersects', 
        target_coordinate_system=coord_system
    )
    if not joined_nuclei.empty:
        # 4. Identify IDs
        # Group A: IDs to REMOVE (infiammazione)
        mask_inf = joined_nuclei['name'].str.contains("infiammazione", case=False, na=False)
        ids_to_remove = joined_nuclei[mask_inf].index.unique()
        # Group B: IDs to ANNOTATE (fibre_trattate)
        mask_treat = joined_nuclei['name'].str.contains("fibre_trattate", case=False, na=False)
        ids_to_treat = joined_nuclei[mask_treat].index.unique()
        # Prevent annotating nuclei that are about to be removed
        ids_to_treat = ids_to_treat.difference(ids_to_remove)
        # 5. Update the Table
        # First, apply annotations to the table (before filtering rows)
        if len(ids_to_treat) > 0:
            is_treated = table.obs[instance_key].isin(ids_to_treat)
            table.obs.loc[is_treated, 'in_treatment'] = True
        # Second, filter the table rows (remove inflammation)
        if len(ids_to_remove) > 0:
            keep_mask = ~table.obs[instance_key].isin(ids_to_remove)
            table = table[keep_mask].copy()
        sdata[table_name] = table
        # 6. Sync Shapes to Table
        # match_element_to_table returns a tuple: ({element_name: element}, table)
        filtered_elements, _ = sd.match_element_to_table(sdata, region_name, table_name)
        # Extract the specific element (GeoDataFrame) from the dictionary and assign it
        sdata[region_name] = filtered_elements[region_name]
        print(f"Processed {region_name}: Removed {len(ids_to_remove)} nuclei. Annotated {len(ids_to_treat)} nuclei.")
    else:
        print(f"Processed {region_name}: No overlap with polygons.")

# check plot ----------------------------------
#
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata3.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('filtered_nuclei', outline=False, outline_alpha=1, outline_width=1, fill_alpha=1
).pl.show(ax = ax, coordinate_systems=blocco, 
    save = f'output_python/roba_da_buttare/{sample_keys[0]}_noinflamed.png'
)
#
sdata2['nuclei_counts'].obs['in_treatment'] = (
    sdata2['nuclei_counts'].obs['in_treatment']
    .astype(str)
    .astype('category')
)

plt.figure(figsize=(50, 50))
ax = plt.gca()

(
    sdata2.pl.render_images('full_image', scale='scale3')
    .pl.render_shapes(
        'filtered_nuclei',
        color='in_treatment',
        outline=False,
        outline_alpha=1,
        outline_width=1,
        fill_alpha=1
    )
    .pl.show(
        ax=ax,
        coordinate_systems=blocco,
        save = f'output_python/roba_da_buttare/{sample_keys[1]}_intreatment.png'
    )
)

#
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata1.pl.render_images('full_image', scale = 'scale3'
).pl.render_shapes('filtered_nuclei', outline=False, outline_alpha=1, outline_width=1, fill_alpha=1
).pl.show(ax = ax, coordinate_systems=blocco, 
    save = f'output_python/roba_da_buttare/{sample_keys[1]}_noinflamed.png'
)

#
# Ensure the column is categorical so the palette maps correctly
sdata1['nuclei_counts'].obs['in_treatment'] = (
    sdata1['nuclei_counts'].obs['in_treatment']
    .astype(str)
    .astype('category')
)

plt.figure(figsize=(50, 50))
ax = plt.gca()

(
    sdata1.pl.render_images('full_image', scale='scale3')
    .pl.render_shapes(
        'filtered_nuclei',
        color='in_treatment',
        outline=False,
        outline_alpha=1,
        outline_width=1,
        fill_alpha=1
    )
    .pl.show(
        ax=ax,
        coordinate_systems=blocco,
        save=f'output_python/roba_da_buttare/{sample_keys[2]}_intreatment.png'
    )
)

# BLOCCO 4 NON HA FUNZIONATO PER IL CAMPIONE C26SMAD23

# Assuming gfp_poly is your GeoDataFrame
# Plot the GeoDataFrame, coloring by the 'name' column
fig, ax = plt.subplots(figsize=(10, 10))
gfp_co.plot(column='name', ax=ax, legend=True, cmap='Set2')

ax.set_title("GFP Polygons by Name")
plt.axis('off') # Turn off axis numbers/frame if desired

# Save the plot
output_filename = "/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/roba_da_buttare/blocco2_gfp_girato_plot.png"
plt.savefig(output_filename, dpi=300, bbox_inches='tight')

print(f"Plot saved to {output_filename}")


# it's working... i need to save all the modifications
for sdata in sdata_list:
    sdata.delete_element_from_disk('filtered_nuclei')
    sdata.write_element('filtered_nuclei')
    sdata.delete_element_from_disk('nuclei_counts')
    sdata.write_element('nuclei_counts')
    sdata.write_element('GFP_poly')


