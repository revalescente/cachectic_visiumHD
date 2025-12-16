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

# To convert from ome.tif to tif with bftools
# # we convert only the first series (aka the full res serie and not the full pyramid)
# and use big tiff to prevent the 4gb limit
#
# ./bftools/bfconvert -series 0 -bigtiff \
#   Fluo_images/overlayed_ome_tif/blocco9.ome.tif \
#   Fluo_images/warped_tif/blocco9.tif
#

blocco = 'blocco2'

img = io.imread(f"/mnt/europa/valerio/Fluo_images/warped_tif/{blocco}_rotated.tif")
# print(img.shape, img.dtype)

# 2nd channel correct, the other two no, fuck.

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

# check on the image values because something fishy
for c in range(img_rescaled.shape[-1]):
    ch_min = img_rescaled[..., c].min()
    ch_max = img_rescaled[..., c].max()
    print(f"Channel {c}: min={ch_min}, max={ch_max}")

# this creates 4 downscaled images as well
img_parsed = Image2DModel.parse(data=img_rescaled, 
            scale_factors=(2, 2, 2), 
            transformations={blocco: Identity()},
            dims=("y", "x", "c")
)  

temp_sd = sd.SpatialData(images={"fluo_image": img_parsed})


samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]
# data

shape_key = 'filtered_nuclei'
table_key = 'nuclei_counts'
sdata_path = '/mnt/europa/valerio/data/zarr_store/old_samples'
sdata_list = []
for sample, data in block_samples.items():
    # subsetting the temp_sd for the single sample
    min_coord = data['min_coordinate']
    max_coord = data['max_coordinate']
    sd1 = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=blocco
    )
    # reading the sdata 
    sdata = sd.read_zarr(f"{sdata_path}/{blocco}_{sample}_storto.zarr")
    # adding the fluo image, cropped as the HE image
    sdata['fluo_image'] = sd1['fluo_image']
    # channels aggregation 
    channel_aggregation = sopa.aggregation.aggregate_channels(sdata, image_key='fluo_image', shapes_key=shape_key, expand_radius_ratio=0, mode='max', no_overlap=False) 
    max_values_vector = channel_aggregation.max(axis=1)
    sdata[table_key].obs['GFP_value'] = max_values_vector
    sdata_list.append(sdata)

zero_cmap = set_zero_in_cmap_to_transparent(cmap="Wistia")
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata_list[2].pl.render_shapes('filtered_nuclei', color = 'GFP_value', cmap = 'Wistia', table_name = 'nuclei_counts'
).pl.show(ax = ax, coordinate_systems=blocco, save = f'output_python/GFP_nuclei_{blocco}.png')


for sdata in sdata_list:
    # saving the objects 
    sdata.write_element('fluo_image')
    sdata.delete_element_from_disk(table_key) 
    sdata.write_element(table_key)

# 
# 
# sample = 'c26foxO'
# min_coord = samples_dict[blocco][sample]['min_coordinate']
# max_coord = samples_dict[blocco][sample]['max_coordinate']
# 
# sd1 = temp_sd.query.bounding_box(
#                 axes=['x', 'y'],
#                 min_coordinate=min_coord,
#                 max_coordinate=max_coord,
#                 target_coordinate_system=blocco
#             )
# 
# 
# # plt.figure(figsize=(20, 20))
# # ax = plt.gca()
# # sd1.pl.render_images('fluo_image', cmap = 'grey', scale = 'scale3').pl.show(ax = ax, coordinate_systems=blocco, save = 'output_python/roba_da_buttare/align_prova2.png')
# 
# # aggiungiamo a sdata e cerchiamo di rasterizzare l'immagine in modo che diventi una label con aggregazione del canale di interesse
# # all'interno dei nuclei. valore massimo tra tutti i channel
# # 
# sdata = sd.read_zarr(f"/mnt/europa/shared/sandri_visiumHD_data/bins/version_1.0.0/{blocco}_{sample}")
# sdata['fluo_image'] = sd1['fluo_image']
# 
# # retry with sopa... easier?
# channel_aggregation = sopa.aggregation.aggregate_channels(sdata, image_key='fluo_image', shapes_key='intissue_008um', expand_radius_ratio=0, mode='max', no_overlap=False)
# max_values_vector = channel_aggregation.max(axis=1)
# 
# sdata['filtered'].obs['GFP_value'] = max_values_vector
# 
# zero_cmap = set_zero_in_cmap_to_transparent(cmap="Wistia")
# 
# plt.figure(figsize=(20, 20))
# ax = plt.gca()
# sdata.pl.render_images('blocco1_hires_image'
# ).pl.render_shapes('intissue_008um', color = 'GFP_value', cmap = zero_cmap, table_name = 'filtered'
# ).pl.show(ax = ax, coordinate_systems=blocco, save = 'output_python/GFP_bins_b1.png')
# 
# # # save all the necessary things
# # 
# # # It's not possible to overwrite, so I need delete from disk then rewrite the element which remain in-memory storage
# # sdata.write_element('fluo_image')
# # sdata.delete_element_from_disk('nuclei_counts') 
# # sdata.write_element('nuclei_counts')
# # 
