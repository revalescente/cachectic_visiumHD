import numpy as np

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
