import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Identity, get_transformation, remove_transformation
import sopa
import spatialdata_plot
import matplotlib.pyplot as plt
from skimage import io
import numpy as np
from py_scripts.utils.utils_fun import read_from_json
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent

# let's display the areas where no expression is detected as transparent
new_cmap = set_zero_in_cmap_to_transparent(cmap="Blues")
new_cmap



img = io.imread("/mnt/europa/valerio/Fluo_images/warped_tif/blocco3.tif")
print(img.shape, img.dtype)

img_rescaled = np.empty_like(img)
for c in range(img.shape[-1]):
    ch = img[..., c]
    ch_min = ch.min()
    ch_max = ch.max()
    # Avoid division by zero if channel is flat
    if ch_max > ch_min:
        img_rescaled[..., c] = 255 * (ch - ch_min) / (ch_max - ch_min)
    else:
        img_rescaled[..., c] = 0  # or ch_min, or 255, as appropriate

# only first channel
img_ch1 = img_rescaled[:, :, 0:1]
# adding to spatialdata

# this creates 4 downscaled images as well
img_parsed = Image2DModel.parse(data=img_rescaled, 
            scale_factors=(2, 2, 2), 
            transformations={'blocco3': Identity()},
            dims=("y", "x", "c")
)  


temp_sd = sd.SpatialData(images={"fluo_image": img_parsed})

plt.figure(figsize=(20, 20))
ax = plt.gca()
temp_sd.pl.render_images('fluo_image', cmap = 'grey', scale = 'scale3').pl.show(ax = ax, coordinate_systems="blocco3", save = 'output_python/aligned_b3.png')


samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

min_coord = samples_dict['blocco2']['c26SMAD23']['min_coordinate']
max_coord = samples_dict['blocco2']['c26SMAD23']['max_coordinate']

# BLOCCO 2 rotto? blocco1 funziona senza problemi il ritaglio. invece blocco2 no, solo che poi inserendolo nel sdata con l'immagine
# full res non si allineano... di conseguenza anche i nuclei? devo verificare la pipeline?? che due coglioni, provo con immagine
# del blocco3 per capire se ho fatto bene una volta oppure se ho sbagliato solo una volta... bo

# per qualche motivo lo zero della y anziche partire da sotto parte da sopra quindi devo invertire i valori della y
# >>> samples_dict['blocco2']['c26']
# {'min_coordinate': [0, 0], 'max_coordinate': [21414, 7200], 'sample_key': 'blocco2_c26'}
# y = (0, 7200) -> y = (20069, 20069-7200) = (20069, 12869), x = (0,21414)
min_coords = [0,12869]
max_coords = [21414, 20069]
# >>> samples_dict['blocco2']['c26murf1']
# {'min_coordinate': [0, 6500], 'max_coordinate': [21414, 12000], 'sample_key': 'blocco2_c26murf1'}
# y = (6500, 12000) -> y = (20069-6500, 20069-12000)
min_coords = [0, 8069]
max_coords = [21414, 13569]
# >>> samples_dict['blocco2']['c26SMAD23']
# {'min_coordinate': [0, 11500], 'max_coordinate': [21414, 20069], 'sample_key': 'blocco2_c26SMAD23'}
# y = (11500, 20069) -> y = (20069-11500, 0)
min_coords = [0, 0]
max_coords = [21414, 8569]


sd1 = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system='blocco2'
            )
            
sd1 = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=[0,0],
                max_coordinate=[21414,10000],
                target_coordinate_system='blocco2'
            )
        

plt.figure(figsize=(20, 20))
ax = plt.gca()
sd1.pl.render_images('fluo_image', cmap = 'grey', scale = 'scale2').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/aligned_b2_c26smad23_v2.png')


sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco2_c26SMAD23.zarr")
sdata['fluo_image'] = sd1['fluo_image']

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images('full_image', scale = 'scale2').pl.render_images('fluo_image', cmap = 'grey', scale = 'scale2').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/aligned_b2_c26smad23_fullefluo.png')

# il problema e' brutto. l'immagine non sembra solo girata ma anche specchiata.

channel_aggregation = sopa.aggregation.aggregate_channels(sdata, image_key='fluo_image', shapes_key='filtered_nuclei', expand_radius_ratio=0, mode='max', no_overlap=False) 
max_values = np.max(channel_aggregation, axis=1)
sdata['nuclei_counts'].obs['GFP_value'] = max_values


fig, ax = plt.subplots(figsize=(8, 6))
# Plot the histogram of the 'max_channel_value' column
ax.hist(max_values, bins=50, color='skyblue', edgecolor='black')
file_path = '/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/max_channel_value_histogram.png'
plt.savefig(file_path, dpi=300, bbox_inches='tight')

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_shapes('filtered_nuclei', color = 'GFP_value', cmap = "Wistia", table_name = 'nuclei_counts').pl.show(ax = ax, coordinate_systems="blocco2", save = 'output_python/GFP_nuclei_b2_c26SMAD23.png')
