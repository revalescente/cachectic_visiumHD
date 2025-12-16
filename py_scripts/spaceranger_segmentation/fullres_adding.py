import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Identity, get_transformation, remove_transformation, Scale, set_transformation
import sopa
import spatialdata_plot
import matplotlib.pyplot as plt
from skimage import io
import numpy as np
from py_scripts.utils.utils_fun import read_from_json
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent


blocco = 'blocco2'

img = io.imread(f"/mnt/europa/valerio/HE_images/color_corrected/pp_{blocco}_20x.tif")

# this creates 4 downscaled images as well
img_parsed = Image2DModel.parse(data=img, 
            scale_factors=(2, 2, 2), 
            transformations={blocco: Identity()},
            dims=("y", "x", "c")
)  

temp_sd = sd.SpatialData(images={"full_image": img_parsed})

samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
block_samples = samples_dict.get(blocco)
sample_keys = [data['sample_key'] for data in block_samples.values()]

# scale factor for the blocco
scalejs = read_from_json(f'/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell/{blocco}/outs/segmented_outputs/spatial/scalefactors_json.json')
scale_factor = 1 / scalejs['tissue_hires_scalef']

# single sample
sdata = sd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples/{sample_keys[0]}")
min_coord = [0, 0]
max_coord = [16166, 8000]
sd1 = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=blocco
)
# adding the fullres image, cropped
sdata['full_image'] = sd1['full_image']

# 3. Create the Scale transformation object
# We apply the scale to both X and Y axes
scaler_transform = Scale([scale_factor, scale_factor], axes=("x", "y"))

# 4. Apply the transformation to the nuclei boundaries
# We map the element to the target coordinate system 'blocco1'
shape_key = 'blocco1_c26STAT3_nuclei_boundaries'

set_transformation(
    element=sdata.shapes[shape_key],
    transformation=scaler_transform,
    to_coordinate_system=blocco
)
# 
# 
plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.pl.render_images(
  'blocco2_hires_tissue_image'
).pl.show(ax = ax, coordinate_systems='downscale_to_hires', save = f'output_python/sr_b2_crying.png')



# all blocco
table_key = 'segmentation_counts'
sdata_path = '/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/samples'
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
    sdata = sd.read_zarr(f"{sdata_path}/{blocco}_{sample}")
    shape_key = f'{blocco}_{sample}_nuclei_boundaries'
    # adding the fullres image, cropped
    sdata['full_image'] = sd1['full_image']
    # define the transformation for the nuclei_boundaries
    scaler_transform = Scale([scale_factor, scale_factor], axes=("x", "y"))
    set_transformation(
        element=sdata.shapes[shape_key],
        transformation=scaler_transform,
        to_coordinate_system=blocco
    )
    # add to sdata_list
    sdata_list.append(sdata)


plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata_list[0].pl.render_images(
  'full_image', scale = 'scale3'
).pl.render_shapes('blocco2_c26_nuclei_boundaries').pl.show(ax = ax, coordinate_systems=blocco, save = f'output_python/sr_b2_crying.png')



for sdata in sdata_list:
    # saving the objects 
    sdata.write_element('full_image')
    sdata.write_metadata()
