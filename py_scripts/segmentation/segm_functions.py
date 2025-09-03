import geopandas as gpd
import pandas as pd
from spatialdata.models import ShapesModel

# Function to spatial join bins with respect to the nuclei
def assign_bins_to_nearest_nucleus(bins, nuclei_bounds):
    """
    For each bin that intersects one or more nuclei, assign it to the nearest nucleus boundary
    (by centroid-to-boundary distance), and return a GeoDataFrame with bin-to-nucleus assignment.
    """
    # Ensure location_id is categorical
    bins = bins.copy()
    bins['location_id'] = bins['location_id'].astype('category')

    # Ensure nuclei_bounds has cell_id as a column
    nuclei_bounds = nuclei_bounds.copy()
    nuclei_bounds['cell_id'] = nuclei_bounds.index

    # Step 1: Spatial join to get all intersecting bin-nucleus pairs
    filtered_bins = gpd.sjoin(
        bins,
        nuclei_bounds[['geometry', 'cell_id']],
        how='inner',
        predicate='intersects'
    )

    # Step 2: Merge the nucleus boundary geometry into filtered_bins
    filtered_bins = filtered_bins.merge(
        nuclei_bounds[['cell_id', 'geometry']].rename(columns={'geometry': 'nucleus_geometry'}),
        on='cell_id'
    )

    # Step 3: Compute bin centroid
    filtered_bins['bin_centroid'] = filtered_bins.geometry.centroid

    # Step 4: Compute distance from bin centroid to nucleus boundary
    filtered_bins['distance'] = filtered_bins.apply(
        lambda row: row['bin_centroid'].distance(row['nucleus_geometry']),
        axis=1
    )

    # Step 5: For each bin, keep only the closest nucleus
    nearest_bins = filtered_bins.loc[
        filtered_bins.groupby('location_id', observed=True)['distance'].idxmin()
    ].reset_index(drop=True)

    # Step 6: Select relevant columns for output
    bins_nuclei = nearest_bins[['location_id', 'cell_id', 'name', 'geometry', 'distance']]
    final = ShapesModel.parse(bins_nuclei)
    return final

# ------------------------------------------------------------------------------

def nuclei_filtering(props_df, filters):
    """
    Filter nuclei features DataFrame by user-specified criteria.
    
    Args:
        props_df (pd.DataFrame): DataFrame with nuclei features.
        filters (dict): Filtering criteria, e.g.:
            {
                'area': (min_area, max_area),
                'eccentricity': (min_ecc, max_ecc),
                'solidity': (min_sol, max_sol),
                'extent': (min_ext, max_ext)
            }
            The value for each key can be a tuple (min, max), or just min/max (None for no limit).
    
    Returns:
        pd.DataFrame: Filtered DataFrame.
    """
    df = props_df.copy()
    for col, val in filters.items():
        if isinstance(val, tuple):
            min_val, max_val = val
        else:
            min_val, max_val = val, None
        if min_val is not None:
            df = df[df[col] >= min_val]
        if max_val is not None:
            df = df[df[col] <= max_val]
    return df

# Example usage:
# filters = {'area': (1000, None), 'eccentricity': (None, 0.7), 'solidity': (0.95, None)}
# filtered = nuclei_filtering(props_df, filters)

# ------------------------------------------------------------------------------

# def features_extraction(sdata, nuclei_element_name = "nuclei_boundaries"):
#   '''
#   Function to extract features from the nuclei shapes
#   '''
#   # check: stop if doesn't check: nuclei_element_name must be in the sdata object
#   # Check if attrs are correct:
#   if sopa_attrs_check(sdata) is True: 
#     continue
#   else print("Check sopa's attributes")  # here should stop
#   
#   # blocks key to call other elements
#   blocco_key = re.search(r'(blocco\d+)', sdata.attrs['boundaries_shapes']).group(1)
#   
#   # Let's rasterize the nuclei polygons (trasform in a label image)
#   # extract the extent to rasterize
#   element_extent = sd.get_extent(sdata[nuclei_element_name], coordinate_system=blocco_key, exact=True)
#   # exemple {'x': (1994.433953335766, 15784.49036300504), 'y': (299.5143364244573, 7999.5)}
# 
#   # rasterize shapes
#   sdata[f"{blocco_key}_raster_nuclei"] = sd.rasterize(
#     sdata[nuclei_element_name],
#     ["x", "y"],
#     min_coordinate=[element_extent['x'][0],element_extent['y'][0]],
#     max_coordinate=[element_extent['x'][1],element_extent['y'][1]],
#     target_coordinate_system=blocco_key,
#     target_unit_to_pixels=1,
#   )
#   # transform into integer values
# 
#   # If your array is an xarray with shape (1, H, W), get the first channel
#   label_mask = sdata[f'{blocco_key}_raster_nuclei'].values
#   
#   # If shape is (1, H, W), squeeze to (H, W)
#   if label_mask.ndim == 3 and label_mask.shape[0] == 1:
#       label_mask = label_mask[0]
#   else label_mask
#   
#   # Cast to integer type
#   label_mask = label_mask.astype(np.int32)
# 
#   # add features of the rasterized nuclei boundaries
#   
#   # 2. Compute regionprops
#   props = regionprops_table(label_mask, properties=[
#       'label', 'area', 'eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length'])
#   props_df = pd.DataFrame(props)
#   
#   # Get mapping dictionary from xarray attributes
#   label_to_id = b1_stat3['raster_nuclei'].attrs['label_index_to_category']
#   
#   # Map: create a new column with the nucleus ID
#   props_df['cell_id'] = props_df['label'].map(label_to_id)


# ------------------------------------------------------------------------------

def sopa_attrs_check(sdata):
    """
    Checks if sdata.attrs contains correct SOPA keys and values
    based on the sdata object's name (format: blocco_expcond).
    """
    # Try to get the name from sdata (adjust as needed for your implementation)
    name = getattr(sdata, 'name', None)
    if name is None:
        # Try to infer from zarr path
        path = getattr(sdata, 'zarr_path', getattr(sdata, 'path', None))
        if path is not None:
            import os
            # Get last part of path, remove extension
            name = os.path.splitext(os.path.basename(path))[0]
        else:
            raise ValueError("Cannot determine name of sdata object.")

    # Extract blocco part (assume format blocco_expcond)
    blocco = name.split('_')[0]
    cell_key = f"{blocco}_full_image"
    tissue_key = f"{blocco}_full_image"
    boundaries_key = f"{blocco}_intissue"
    bins_key = "filtered"

    required = {
        "cell_segmentation_key": cell_key,
        "tissue_segmentation_key": tissue_key,
        "bins_table_key": bins_key,
        "boundaries_key": boundaries_key,
    }
    attributes = sdata.attrs
    for key, expected_value in required.items():
        if key not in attributes:
            print(f"Missing attribute: {key}")
            return False
        if attributes[key] != expected_value:
            print(f"Attribute {key} has value {attributes[key]}, expected {expected_value}")
            return False
    return True

# ------------------------------------------------------------------------------

# def apply_segmentation(sdata, patch_width = 500, patch_overlap = 50, stardist_model = '2D_veratile_he',
# nuclei_key = "nuclei_boundaries", **kwargs):
#   '''
# 
#   Function to apply segmentation to out data, one sample at a time
# 
#   attrs structure:
#     sdata.attrs:
#       {'bins_table': 'filtered',
#        'boundaries_shapes': 'blocco9_intissue',
#        'cell_segmentation_image': 'blocco9_full_image',
#        'tissue_segmentation_image': 'blocco9_full_image'}
# 
# 
#   '''
#   # Check if attrs are correct:
#   if sopa_attrs_check(sdata) is True:
#     continue
#   else print("Check sopa's attributes")  # here should stop
# 
#   # key to match sdata elements
#   blocco_key = re.search(r'(blocco\d+)', sdata.attrs['boundaries_shapes']).group(1)
#   nuclei_element_name = f'{blocco_key}_{nuclei_key}'
# 
#   # create overlapping patches
#   sopa.make_image_patches(sdata, patch_width = patch_width, patch_overlap = patch_overlap, **kwargs)
# 
#   # segmentation
#   sopa.segmentation.stardist(sdata, model_type=stardist_model,
#   image_key= sdata.attrs['cell_segmentation_image'], min_area=10, delete_cache=True,
#   recover=False, prob_thresh=0.2, nms_thresh=0.6, key_added = nuclei_element_name, **kwargs)
# 
#   sdata[f'{blocco_key}_bins_nuclei'] = assign_bins_to_nearest_nucleus(bins = sdata[f'{blocco_key}_intissue_002um'],
#   nuclei_bounds = sdata[nuclei_element_name])
# 
#   ... to continue ...

  



































