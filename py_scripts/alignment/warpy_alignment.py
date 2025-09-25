# PyImageJ / Scyjava
from scyjava import jimport
import imagej
import numpy as np
# JPype
from jpype.types import JString, JArray
from skimage import io
import matplotlib.pyplot as plt
import py_scripts.alignment.align_funs as af
# from importlib import reload
# reload(sf)


def get_java_dependencies():
    """
    Returns the jar files that need to be included into the classpath
    :return:
    """
    return [# 'net.imagej:imagej:2.16.0',
           	'ch.epfl.biop:bigdataviewer-biop-tools:0.11.2'
            # add another jar here if necessary
    ]
# Start a Fiji instance - it's big and will take time the first time.
ij = imagej.init(get_java_dependencies(), mode="headless")
# Warpy uses a json serialization of its real transform objects. They can be accessed and used in python through pyimagej
transform_file_path = '/mnt/europa/valerio/repositories/cachetic_visiumHD/json/transform_json/transform_10_8_b1_c26foxO.json'
# Importing Java classes in python
RealTransformHelper = jimport('net.imglib2.realtransform.RealTransformHelper')
# Retrieve the transformation object
transform = RealTransformHelper.fromJson(ij.context(), JString(transform_file_path))
# Transforming the fluo image by using RealTransformHelper ----------

# Read a TIFF image
img = io.imread("/mnt/europa/valerio/Fluo_images/samples/blocco1_c26foxO.tif")
# img2 = io.imread("/mnt/europa/valerio/HE_images/color_corrected/pp_blocco1_20x.tif")

# Transform the fluo image.
RealTransformHelper.apply(ij.py.to_java(img), transform)

# save it
io.imsave("/mnt/europa/valerio/Fluo_images/warped/blocco1_c26foxO.tif", img)

plt.figure(figsize = (20,20))
plt.imshow(img)  # Use cmap='gray' for grayscale images
plt.axis('off')  # Hide axes
plt.savefig('/mnt/europa/valerio/figures/warpy_b1c26foxO.png', bbox_inches='tight', pad_inches=0)
plt.close()

# plot of the images superimposed
fig, ax = plt.subplots(figsize=(20, 20))
# Plot first image
ax.imshow(img2, alpha=0.5)  # Change alpha as needed
# Plot second image on top
ax.imshow(img)  # Use alpha to blend images
ax.axis('on')
plt.savefig('/mnt/europa/valerio/figures/align_1.png', bbox_inches='tight', pad_inches=0)
plt.close()

# ------------------------------------------------------------------------------
# With functions

img = io.imread("/mnt/europa/valerio/Fluo_images/samples/blocco1_c26foxO.tif")

modded = af.transform_to_physical_coordinates(img, pixel_size_x = 0.315, pixel_size_y = 0.315)





















