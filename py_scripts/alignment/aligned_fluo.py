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

# To convert from ome.tif to tif with bftools
# # we convert only the first series (aka the full res serie and not the full pyramid)
# and use big tiff to prevent the 4gb limit
#
# ./bftools/bfconvert -series 0 -bigtiff \
#   Fluo_images/overlayed_ome_tif/blocco9.ome.tif \
#   Fluo_images/warped_tif/blocco9.tif
#

img = io.imread("/mnt/europa/valerio/Fluo_images/blocco1.tif")
print(img.shape, img.dtype)

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
  
# second channel only

img_ch2 = img_rescaled[:, :,1:2]

# this creates 4 downscaled images as well
img_parsed = Image2DModel.parse(data=img_ch2, 
            scale_factors=(2, 2, 2), 
            transformations={'blocco1': Identity()},
            dims=("y", "x", "c")
)  

# spe_b1 = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/filtered/filtered_blocco1.zarr")
# spe_b1['fluo_image'] = img_ch2

plt.figure(figsize=(20, 20))
ax = plt.gca()
spe_b1.pl.render_images('fluo_image2', ["green","green","green"], scale = "scale2").pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/aligned_maybe3.png')

#' ora ho l'immagine in un sdata, devo inserirla all'interno dei rispettivi campioni,
#' credo di dover definire la trasformazione spaziale identica all'immagine full_image
#' che uso nei vari campioni, così dovrei riuscire ad inserirla già allineata.
#' 

# trans = get_transformation_from_(sdata['full_image'], 'blocco1')
# remove_transformation(img_parsed, 'blocco1')
# 
# set.transformation
# 
# img_parsed = Image2DModel.parse(data=img_parsed, 
#             scale_factors=(2, 2, 2),
#             transformations = trans,
#             dims=("y", "x", "c")
# )  

# facciamo cosi, definisco un oggetto sd con solo l'immagine, la taglio come ho tagliato gli altri campioni e poi copio l'immagine risultante
# all'interno del rispettivo campione. 

temp_sd = sd.SpatialData(images={"fluo_image": img_parsed})

samples_dict = read_from_json('/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json')

min_coord = samples_dict['blocco1']['c26foxO']['min_coordinate']
max_coord = samples_dict['blocco1']['c26foxO']['max_coordinate']

sd1 = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system='blocco1'
            )


plt.figure(figsize=(20, 20))
ax = plt.gca()
sd1.pl.render_images('fluo_image', cmap = 'grey').pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/aligned_sub_b1_v2.png')

# aggiungiamo a sdata e cerchiamo di rasterizzare l'immagine in modo che diventi una label con aggregazione del canale di interesse
# all'interno dei nuclei. proviamo intanto con un valore medio del channel corrispondente

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco1_c26foxO.zarr")
sdata['fluo_image'] = sd1['fluo_image']

#' spatialdata.aggregation: adesso aggrego i valori dell'immagine a fluorescenza rispetto ai nuclei segmentati. quindi
#' 
#' 
#' rasterizzo nuclei
# 
# rasterized_raw_image = rasterize(
#     sdata["raw_image"],
#     min_coordinate=[-2000, -2000],
#     max_coordinate=[10_000, 10_000],
#     axes=("x", "y"),
#     target_coordinate_system="global",
#     target_unit_to_pixels=1,
# )
# sdata["rasterized_raw_image"] = rasterized_raw_image
# 
# sdata.pl.render_images("rasterized_raw_image").pl.render_shapes(fill_alpha=0.2).pl.show()
# 
# sd.get_extent(sdata['filtered_nuclei'], coordinate_system = 'blocco1')
# {'x': (2568.5152483790284, 14136.484657115801), 'y': (15370.503506204665, 21060.4846571158)}
# 
# rast_nuclei = sd.rasterize('filtered_nuclei', axes = ('x', 'y'), 
#                 min_coordinate = [2568.5152483790284, 15370.503506204665], 
#                 max_coordinate = [14136.484657115801, 21060.4846571158], 
#                 target_coordinate_system = 'blocco1', target_unit_to_pixels=1, sdata=sdata
# )
# sdata.labels['rast_nuclei'] = rast_nuclei
# rast_nuclei_pars = Labels2DModel.parse(data=rast_nuclei, 
#             dims=("y", "x")
# )  
# sdata_im = sdata.aggregate(values="fluo_image", by="rast_nuclei", agg_func="mean", target_coordinate_system = 'blocco1')
# 

# retry with sopa... easier?
channel_aggregation = sopa.aggregation.aggregate_channels(sdata, image_key='fluo_image', shapes_key='filtered_nuclei', expand_radius_ratio=0, mode='max', no_overlap=False) 

sdata['nuclei_counts'].obs['GFP_value'] = channel_aggregation

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_shapes('filtered_nuclei', color = 'GFP_value', cmap = "Wistia", table_name = 'nuclei_counts').pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/GFP_nuclei.png')


plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=[8500, 16000],
                max_coordinate=[11000, 18000],
                target_coordinate_system='blocco1'
).pl.render_shapes('filtered_nuclei', color = 'GFP_value', cmap = "Wistia", table_name = 'nuclei_counts').pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/GFP_nuclei_zoom.png')

# funziona in effetti.

# save all the necessary things

# It's not possible to overwrite, so I need delete from disk then rewrite the element which remain in-memory storage
sdata.write_element('fluo_image')
sdata.delete_element_from_disk('nuclei_counts') 
sdata.write_element('nuclei_counts')

# ------------------------------------------------------------------------------

# filtro piu nuclei perche ce ne sono troppi troppo grandi

areas = sdata['nuclei_counts'].obs['area']

plt.figure(figsize=(10, 6))
plt.hist(areas, bins=50, color='steelblue', edgecolor='black', alpha=0.7)

# Add statistics lines
plt.axvline(areas.mean(), color='red', linestyle='dashed', linewidth=1, 
           label=f'Mean: {areas.mean():.2f}')
plt.axvline(areas.median(), color='green', linestyle='dashed', linewidth=1, 
           label=f'Median: {areas.median():.2f}')

# Add labels and title
plt.xlabel('Nucleus Area', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.title('Distribution of Nucleus Areas', fontsize=14)
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

# Save the figure
plt.savefig('/mnt/europa/valerio/repositories/cachetic_visiumHD/figures/output_python/nucleus_area_histogram.png', dpi=300)  # PNG format

# -------------------------------------------------------------------------

# filter nuclei too big... but how much??? I'll follow the path of spaceranger

# Filter based on a numerical threshold
adata_filt = sdata['nuclei_counts'][sdata['nuclei_counts'].obs['area']<=1500]
sdata['filter_try'] = adata_filt

plt.figure(figsize=(50, 50))
ax = plt.gca()
sdata.pl.render_shapes('filtered_nuclei', color = 'GFP_value', cmap = "Wistia", table_name = 'filter_try').pl.show(ax = ax, coordinate_systems="blocco1", save = 'output_python/GFP_nuclei_filtered.png')
