import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
from skimage import io
import geopandas as gpd
import pandas as pd
from spatialdata.models import ShapesModel, Image2DModel
from spatialdata.transformations import Identity, set_transformation, get_transformation
from py_scripts.utils.utils_fun import read_from_json

samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")

blocco = 'blocco1'
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]
# sopa
# sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[0]}.zarr")
# sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[1]}.zarr")
# sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/samples/{sample_keys[2]}.zarr")
# spaceranger
sdata1 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[0]}")
sdata2 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[1]}")
sdata3 = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[2]}")
sdata_list = [sdata1, sdata2, sdata3]


# all polygons
geojson_path = f"/mnt/europa/valerio/data/json/geojson_dir/in_tissue/fullres/{blocco}_allpolys_fullres.geojson" 


poly = gpd.read_file(geojson_path)
poly = poly.set_crs(None, allow_override=True)
# proviamo a separare i polygoni cosi matplot fa meno fatica... ???
#gfp_poly = gfp_poly.explode(index_parts=False).reset_index(drop=True)
poly_parse = ShapesModel.parse(poly, transformations={blocco: Identity()})

sdata1.shapes["all_poly"] = poly_parse
sdata2.shapes["all_poly"] = poly_parse
sdata3.shapes["all_poly"] = poly_parse

# combined plot
# Initialize the figure and axis
plt.figure(figsize=(30, 30))
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
        'all_poly', 
        outline=False, 
        outline_alpha=1, 
        outline_width=1, 
        fill_alpha=0
    ).pl.show(
        ax=ax, 
        coordinate_systems=blocco  # Argument for the coordinate system
    )

# Save the combined plot
plt.savefig(f'/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/roba_da_buttare/{blocco}_combined_poly.png')
plt.show()
