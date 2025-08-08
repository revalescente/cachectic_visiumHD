import os
import numpy as np
from astropy.io import fits
from skimage import io, color
from scipy.signal import convolve2d
from tqdm import tqdm

# HSV color conversion and keep saturation only
def extract_saturation_channel(rgb_image):
    """
    Extract saturation channel from HSV colors
    """
    hsv_image = color.rgb2hsv(rgb_image)
    return hsv_image[..., 1]  # Canale S (saturation)
  
# invert and sqrt one channel
def process_saturation_channel(saturation):
    """
    apply inversion and sqrt transformation
    """
    inverted = 1.0 - saturation  # Inversione
    sqrt = np.sqrt(inverted)     # Radice di 2
    return sqrt


def load_kernel(kernel_path):
    """Carica il kernel FITS e lo normalizza"""
    hdul = fits.open(kernel_path)
    kernel = hdul[0].data.astype(np.float32)
    hdul.close()
    return kernel / np.sum(kernel)  # Normalizza per conservare la luminosità

def process_convolutions(input_dir, output_dir, kernel_path):
    """
    Applica la convoluzione a tutte le immagini nella cartella input_dir
    usando il kernel specificato e salva i risultati in output_dir
    """
    # Carica il kernel
    kernel = load_kernel(kernel_path)
    
    # Crea la cartella di output se non esiste
    os.makedirs(output_dir, exist_ok=True)
    
    # List of images to process (with prefix 'pp_' from the visual_step.py)
    image_files = [f for f in os.listdir(input_dir) 
                  if f.startswith('pp_') and f.lower().endswith(('.tif', '.tiff'))]
    
    for filename in tqdm(image_files, desc="Images convolution"):
        # load image
        image_path = os.path.join(input_dir, filename)
        rgb_image = io.imread(image_path)
        
        # convert in float [0,1] if it isn't uint8
        if rgb_image.dtype == np.uint8:
            rgb_image = rgb_image.astype(np.float32) / 255.0
        
        # Apply the tranformations
        saturation = extract_saturation_channel(rgb_image)
        processed = process_saturation_channel(saturation)
        
        # Apply the convolution
        convolved = convolve2d(processed, kernel, mode='same', boundary='symm')
        
        # save the output
        output_filename = f"convolved_{filename.replace('pp_', '')}"
        output_path = os.path.join(output_dir, output_filename)
        
        # 8-bit normalization
        convolved_normalized = ((convolved - np.min(convolved)) / (np.max(convolved) - np.min(convolved)) * 255).astype(np.uint8)
        io.imsave(output_path, convolved_normalized)

# --- Configuration ---
if __name__ == "__main__":
    # Folder with the color corrected images
    input_dir = "HE_images/color_corrected"
    
    # Output folder
    output_dir = "HE_images/preprocessed"
    
    # kernel FITS folder
    kernel_path = "data/sinc_kernel/sinc1_fits.fit"
    
    # process the images
    process_convolutions(input_dir, output_dir, kernel_path)
    print("Work done!")
