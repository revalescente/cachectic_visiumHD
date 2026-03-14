import zarr
import os
import shutil

# PATH TO YOUR ZARR STORE
ZARR_PATH = "/mnt/europa/valerio/data/zarr_store/samples/blocco2_c26SMAD23.zarr"

def check_group(group_path, element_type):
    """
    Checks a specific group (images or labels) for valid OME-NGFF metadata.
    """
    if not os.path.exists(group_path):
        return
    print(f"\nScanning {element_type} in {group_path}...")
    # Open the group using standard Zarr (safe mode)
    try:
        grp = zarr.open(group_path, mode='r')
    except Exception as e:
        print(f"  ⚠️  Could not open group: {e}")
        return
    # Iterate over all elements (images/labels) inside
    for key in grp.keys():
        element_path = os.path.join(group_path, key)
        sub_grp = grp[key]
        
        status = "✅ OK"
        is_broken = False
        # Check 1: Does it have attributes?
        if not hasattr(sub_grp, 'attrs'):
             status = "❌ BROKEN (No attributes)"
             is_broken = True
        # Check 2: Does it have 'multiscales'?
        elif 'multiscales' not in sub_grp.attrs:
             status = "❌ BROKEN (Missing 'multiscales' metadata)"
             is_broken = True
        print(f"  - {key}: {status}")
        if is_broken:
            print(f"    COMMAND TO FIX:  rm -rf {element_path}")
# Run the check
if __name__ == "__main__":
    if not os.path.exists(ZARR_PATH):
        print("Zarr path does not exist.")
        exit(1)
    print(f"Inspecting: {ZARR_PATH}")
    check_group(os.path.join(ZARR_PATH, "images"), "IMAGES")
    check_group(os.path.join(ZARR_PATH, "labels"), "LABELS")
    check_group(os.path.join(ZARR_PATH, "shapes"), "SHAPES") # Sometimes shapes cause issues too
