import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Identity, get_transformation, remove_transformation
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent
import spatialdata_plot
import matplotlib.pyplot as plt
import sopa
from skimage import io
import numpy as np
import os
from py_scripts.utils.utils_fun import read_from_json

# let's redo it with blocco 2
img = io.imread("/mnt/europa/valerio/Fluo_images/warped_tif/blocco2_try.tif")

print(img.shape, img.dtype)

# check on the image values because something fishy
num_channels = img.shape[2]
for i in range(num_channels):
    # Select the current channel
    channel_data = img[:, :, i]
    # Calculate min and max
    min_val = channel_data.min()
    max_val = channel_data.max()
    print(f"Channel {i}:  Min = {min_val},  Max = {max_val}")
    
# Channel 0:  Min = 0.0,  Max = 255.0
# Channel 1:  Min = 0.0,  Max = 77.0
# Channel 2:  Min = 0.0,  Max = 81.0
    
output_dir = "/mnt/europa/valerio/figures/"  
image_path = "/mnt/europa/valerio/Fluo_images/warped_tif/blocco1_try2.tif"

if img.ndim == 3:
    # Get the base name of the file to use in the output filenames
    base_filename = os.path.splitext(os.path.basename(image_path))[0]
    # Get the number of channels
    num_channels = img.shape[2]
    for i in range(num_channels):
        print(f"Processing Channel {i}...")
        # Isolate the current channel
        channel_data = img[:, :, i]
        # Create a new figure
        plt.figure(figsize=(10, 10))
        # Display the channel using a grayscale colormap
        im = plt.imshow(channel_data, cmap='gray')
        # Add a colorbar to show the intensity scale
        plt.colorbar(im)
        # Add a title
        plt.title(f"{base_filename} - Channel {i}", fontsize=14)
        # Hide the axes for a cleaner look
        plt.axis('off')
        # Define the full path for the output file
        output_filepath = os.path.join(output_dir, f"{base_filename}_channel_{i}.png")
        # Save the figure
        plt.savefig(output_filepath, dpi=150, bbox_inches='tight')
        # Close the figure to free up memory before the next loop
        plt.close()
        print(f"  -> Saved plot to: {output_filepath}")
    print("\nAll channels have been plotted and saved.")
else:
    print("The loaded image is not a 3D multi-channel image. No plots were generated.")


# another try ------------------------------------------------------------------
from bioio import BioImage
import bioio_ome_tiff
import numpy as np
import tifffile
import os
import matplotlib.pyplot as plt

img = BioImage("/mnt/europa/valerio/Fluo_images/overlayed_ome_tif/blocco1.ome.tif", reader=bioio_ome_tiff.Reader)
img.data

channel_data = np.squeeze(img.data)
print(f"New shape (CYX): {channel_data.shape}")

output_path = "/mnt/europa/valerio/Fluo_images/warped_tif/blocco1_from_bioio.tif"

# Use tifffile.imwrite to save the NumPy array.
# It correctly handles the (channels, height, width) shape.
print(f"Saving the 3-channel image to: {output_path}")
tifffile.imwrite(output_path, channel_data)

tifffile.imread(output_path)

output_dir = "/mnt/europa/valerio/figures/"  
if img.ndim == 3:
    # Get the base name of the file to use in the output filenames
    base_filename = os.path.splitext(os.path.basename(output_path))[0]
    # Get the number of channels
    num_channels = img.shape[0]
    for i in range(num_channels):
        print(f"Processing Channel {i}...")
        # Isolate the current channel
        channel_data = img[i, :, :]
        # Create a new figure
        plt.figure(figsize=(10, 10))
        # Display the channel using a grayscale colormap
        im = plt.imshow(channel_data, cmap='gray')
        # Add a colorbar to show the intensity scale
        plt.colorbar(im)
        # Add a title
        plt.title(f"{base_filename} - Channel {i}", fontsize=14)
        # Hide the axes for a cleaner look
        plt.axis('off')
        # Define the full path for the output file
        output_filepath = os.path.join(output_dir, f"{base_filename}_channel_{i}.png")
        # Save the figure
        plt.savefig(output_filepath, dpi=150, bbox_inches='tight')
        # Close the figure to free up memory before the next loop
        plt.close()
        print(f"  -> Saved plot to: {output_filepath}")
    print("\nAll channels have been plotted and saved.")
else:
    print("The loaded image is not a 3D multi-channel image. No plots were generated.")



# ---------------------------------------------------------------------------------

plt.figure(figsize=(10, 10))
# Display the channel using a grayscale colormap
im = plt.imshow(img, cmap='gray')
# Add a colorbar to show the intensity scale
plt.colorbar(im)
# Add a title
# Hide the axes for a cleaner look
plt.axis('off')
# Define the full path for the output file
output_filepath = os.path.join(output_dir, "blocco1_mergedchs.png")
# Save the figure
plt.savefig(output_filepath, dpi=150, bbox_inches='tight')
# Close the figure to free up memory before the next loop
plt.close()

tifffile.imwrite(output_path, img_merged)
    
  

# ----------------------------------------------------------------------------
import os
import numpy as np
from skimage import io, exposure
import matplotlib.pyplot as plt

channel_0_data = img_3ch[0, :, :]

p2, p98 = np.percentile(channel_0_data, (2, 98))
img_enhanced = exposure.rescale_intensity(channel_0_data, in_range=(p2, p98))

print("Contrast enhancement for channel 0 complete.")


# --- 4. Plot and Save the Comparison ---
print("Generating and saving plots...")

# Create a figure with two subplots, side-by-side
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
fig.suptitle("Comparison for Channel 0", fontsize=20)


# Plot 1: The Original Channel 0
ax = axes[0]
im = ax.imshow(channel_0_data, cmap='gray')
ax.set_title('Channel 0 - Original', fontsize=16)
ax.axis('off')
fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)


# Plot 2: The Enhanced Channel 0
ax = axes[1]
im_enhanced = ax.imshow(img_enhanced, cmap='gray')
ax.set_title('Channel 0 - Enhanced Contrast', fontsize=16)
ax.axis('off')
fig.colorbar(im_enhanced, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)


# Adjust layout and save the figure
plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust for suptitle
output_dir = "/mnt/europa/valerio/figures/"
output_filepath = os.path.join(output_dir, "blocco1_channel_0_enhanced.png")
plt.savefig(output_filepath, dpi=200, bbox_inches='tight', facecolor='white')
plt.close() # Close the figure to free memory

print(f"\nSuccessfully saved comparison plot to:\n{output_filepath}")


# --------------------------------------------------------------------------------


try:
    img = io.imread(image_path)
    print(f"Successfully loaded 3-channel image from: {image_path}")
    print(f"Original shape: {img.shape}")
    print(f"Original data type: {img.dtype}\n")
except FileNotFoundError:
    print(f"Error: The file was not found at {image_path}")
    exit()

# --- 2. Isolate Channels ---
if img.ndim != 3 or np.min(img.shape) < 3:
    print("Error: This does not appear to be a 3-channel image.")
    exit()

channel_axis = np.argmin(img.shape)
print(f"Assuming channel axis is dimension {channel_axis}.")

if channel_axis == 0:  # Shape is (C, H, W)
    ch0, ch1, ch2 = img[0, :, :], img[1, :, :], img[2, :, :]
else:  # Shape is (H, W, C)
    ch0, ch1, ch2 = img[:, :, 0], img[:, :, 1], img[:, :, 2]

# --- 3. Normalize Channels 1 and 2 ---
# `exposure.rescale_intensity` scales the image's min/max values to a new range.
# We will scale them to the 0-255 range and convert them to an 8-bit integer type (uint8).

print("Normalizing channels 1 and 2 to the 0-255 range...")

# Normalize Channel 1
ch1_normalized = exposure.rescale_intensity(ch1, out_range=(0, 255)).astype(np.uint8)

# Normalize Channel 2
ch2_normalized = exposure.rescale_intensity(ch2, out_range=(0, 255)).astype(np.uint8)

print("Normalization complete.")
print(f"Channel 0 dtype: {ch0.dtype}, Range: {ch0.min()}-{ch0.max()}")
print(f"Channel 1 new dtype: {ch1_normalized.dtype}, Range: {ch1_normalized.min()}-{ch1_normalized.max()}")
print(f"Channel 2 new dtype: {ch2_normalized.dtype}, Range: {ch2_normalized.min()}-{ch2_normalized.max()}\n")


# --- 4. Re-combine into a New 3-Channel Image for Plotting ---
# Matplotlib's imshow expects the channel dimension to be the last one (H, W, C).
# We use np.stack to join the arrays along a new axis.
print("Stacking channels into a new composite image...")
composite_img = np.stack([ch0, ch1_normalized, ch2_normalized], axis=2)
print(f"New composite image shape: {composite_img.shape}")


# --- 5. Plot and Save the Composite Image ---
print("Generating and saving the composite plot...")
plt.figure(figsize=(15, 15))

# When imshow receives a 3D array with a final dimension of 3,
# it automatically interprets it as an RGB image.
plt.imshow(composite_img)

plt.title('Composite Image (Ch0: Original, Ch1/Ch2: Normalized)', fontsize=16)
plt.axis('off') # Hide axes for a cleaner look

# Save the figure
output_filepath = os.path.join(output_dir, "blocco1_composite_image.png")
plt.savefig(output_filepath, dpi=200, bbox_inches='tight', facecolor='black')
plt.close()

print(f"\nSuccessfully saved composite plot to:\n{output_filepath}")


#-------------------------------------------------------------------------------

#' il primo canale e' corretto tra 0 e 255, gli altri due invece sono con max 77 e 81, 
#' voglio riscalare questi due canali cosi l'RGB risultante funziona ugualmente!!

from skimage import exposure

img_rescaled = img.copy()

print("Rescaling channels 1 and 2...")

# --- Rescale Channel 1 ---
# Select channel 1
ch1 = img[:, :, 1]
# Rescale its intensity from its current range (0-77) to the new range (0-255)
# and convert it to an 8-bit integer, which is standard for the 0-255 range.
ch1_rescaled = exposure.rescale_intensity(ch1, out_range=(0, 255)).astype(np.uint8)


# --- Rescale Channel 2 ---
# Select channel 2
ch2 = img[:, :, 2]
# Rescale its intensity from its current range (0-81) to the new range (0-255)
ch2_rescaled = exposure.rescale_intensity(ch2, out_range=(0, 255)).astype(np.uint8)


# --- Update the new image with the rescaled channels ---
# Note: We ensure channel 0 is also uint8 for consistency
img_rescaled[:, :, 0] = img[:, :, 0].astype(np.uint8) 
img_rescaled[:, :, 1] = ch1_rescaled
img_rescaled[:, :, 2] = ch2_rescaled

# --- Verification ---
# Let's check the ranges of the new image to confirm it worked.
print("\nVerifying ranges of the new rescaled image:")
num_channels_rescaled = img_rescaled.shape[2]
for i in range(num_channels_rescaled):
    channel_data = img_rescaled[:, :, i]
    min_val = channel_data.min()
    max_val = channel_data.max()
    print(f"Channel {i}:  Min = {min_val},  Max = {max_val}")


plt.figure(figsize=(15, 15))

# When imshow receives a 3D array with a final dimension of 3,
# it automatically interprets it as an RGB image.
plt.imshow(img_rescaled)

plt.title('immagine riscalata con i canali strani', fontsize=10)
plt.axis('off') # Hide axes for a cleaner look

# Save the figure
output_dir = "/mnt/europa/valerio/figures/"
output_filepath = os.path.join(output_dir, "blocco1_image_rescaled.png")
plt.savefig(output_filepath, dpi=200, bbox_inches='tight', facecolor='black')
plt.close()
