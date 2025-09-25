import numpy as np

# PyImageJ / Scyjava
from scyjava import jimport
import imagej
import numpy as np

# JPype
from jpype.types import JString, JArray
from skimage import io

def transform_to_physical_coordinates(image, pixel_size_x=1.0, pixel_size_y=1.0, offset_x=0, offset_y=0):
    """
    Convert image array to physical coordinates using pixel size.
    
    Parameters
    ----------
    image : ndarray
        Input image array
    pixel_size_x, pixel_size_y : float
        Physical size of each pixel (e.g., μm per pixel)
    offset_x, offset_y : float
        Coordinate system origin offset
        
    Returns
    -------
    ndarray
        Array of (x, y, value) coordinates in physical units for all pixels
    """
    height, width = image.shape[:2]
    coordinates = []
    
    for y in range(height):
        for x in range(width):
            # Convert to physical coordinates
            physical_x = x * pixel_size_x + offset_x
            physical_y = y * pixel_size_y + offset_y
            
            # Get pixel value
            pixel_value = image[y, x]
            
            # Include all pixels regardless of value
            coordinates.append((physical_x, physical_y, pixel_value))
    
    return np.array(coordinates)

# ------------------------------------------------------------------------------

def get_java_dependencies():
    """
    Returns the jar files that need to be included into the classpath
    :return:
    """
    return [# 'net.imagej:imagej:2.16.0',
           	'ch.epfl.biop:bigdataviewer-biop-tools:0.11.2'
            # add another jar here if necessary
    ]

# ------------------------------------------------------------------------------

def images_alignment(
  json_transform = '/mnt/europa/valerio/repositories/cachetic_visiumHD/json/transform_json/transform_10_8_b1_c26foxO.json',
  image_to_transform = None):
  
  # Start a Fiji instance (download 20mb circa only one time)
  ij = imagej.init(get_java_dependencies(), mode="headless")
  
  # Warpy uses a json serialization of its real transform objects. 
  # They can be accessed and used in python through pyimagej
  transform_file_path = json_transform
  # Importing Java classes in python
  RealTransformHelper = jimport('net.imglib2.realtransform.RealTransformHelper')
  # Retrieve the transformation object
  transform = RealTransformHelper.fromJson(ij.context(), JString(transform_file_path))
  # Transforming the fluo image by using RealTransformHelper ----------
  # Read a TIFF image
  # Transform the fluo image.
  image_transformed = RealTransformHelper.apply(ij.py.to_java(image_to_transform), transform)
  return image_transformed







