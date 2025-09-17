import matplotlib.pyplot as plt
from skimage import io
import json
import os

img_path = "C:/Users/Valerio/nuclei_segmentation_cachetic/Fluo_images/composed/BLOCCO2_RGB.tif"
out_dir = "C:/Users/Valerio/nuclei_segmentation_cachetic/Fluo_images/samples"

# Load your fluo_dict (replace with your actual file path if needed)
with open("C:/Users/Valerio/nuclei_segmentation_cachetic/blocco_sample_bbox_FLUO_dict.json", "r") as f:
    fluo_dic = json.load(f)

# --- User selects blocco and sample ---
print("Available blocchi:", list(fluo_dict.keys()))
blocco = input("Select blocco: ")

print("Available samples in", blocco, ":", list(fluo_dict[blocco].keys()))
sample = input("Select sample: ")

#
img = io.imread(img_path)

# Show full image and ask user to select two y-coords
fig, ax = plt.subplots()
ax.imshow(img)
ax.set_title("Click two points on the left or right edge for y1 and y2")
pts = plt.ginput(2)  # User clicks two points
plt.close(fig)

y_coords = sorted([int(pt[1]) for pt in pts])  # y is the second value (vertical axis)
y_min, y_max = y_coords[0], y_coords[1]

# Crop the image (entire x, from y1 to y2)
cropped_img = img[y_min:y_max, :]

# Show cropped result before saving
fig2, ax2 = plt.subplots()
ax2.imshow(cropped_img)
ax2.set_title(f"Cropped image: y={y_min} to y={y_max}\nClose window to continue")
plt.show()

# --- Set x-coords to full image width ---
x_min = 0
x_max = img.shape[1]

# --- Update dict ---
fluo_dict[blocco][sample]["min_coordinate"] = [x_min, y_min]
fluo_dict[blocco][sample]["max_coordinate"] = [x_max, y_max]

print(f"Updated coordinates for {blocco}_{sample}:")
print("min_coordinate:", fluo_dict[blocco][sample]["min_coordinate"])
print("max_coordinate:", fluo_dict[blocco][sample]["max_coordinate"])

filename = f"{blocco}_{sample}.tif"
out_path = os.path.join(out_dir, filename)
io.imsave(out_path, cropped_img)

# ---------------------------
# save the json file with the info 
with open("C:/Users/Valerio/nuclei_segmentation_cachetic/blocco_sample_bbox_FLUO_dict.json", "w") as f:
    json.dump(fluo_dict, f, indent=2)