import geopandas as gpd
import pandas as pd
from spatialdata.models import ShapesModel
import re
import sopa
import time
import numpy as np
from pathlib import Path
import pandas as pd
import spatialdata as sd
from skimage.measure import regionprops_table
from typing import List
import matplotlib.pyplot as plt
from scipy import sparse
from shapely.geometry import Point

# Function to create the vector with the cell_id of the bins 
def bin_to_cell_id_vector(sdata, table_key: str):
    """
    Assigns each bin (column) to its cell_id using the sparse matrix and cell_id array
    from the given sdata and table_key.
    
    Must have the bins_assignments in the obsm of the table corresponding to the table_key

    Parameters
    ----------
    sdata : SpatialData object
        The SpatialData object containing the tables.
    table_key : str
        Key to access the AnnData table in sdata.

    Returns
    -------
    np.ndarray
        Vector of cell_id strings for each bin.
    """
    bins_assignments = sdata[table_key].obsm["bins_assignments"]
    cell_ids = sdata[table_key].obs["cell_id"]        # cell ids to identify bin assignment
    
    argmax = bins_assignments.argmax(axis=0).A1
    mapped_bin = bins_assignments.getnnz(axis=0)
    
    bin_cell_vector = np.where(mapped_bin, cell_ids.values[argmax], None)
    
    return bin_cell_vector

# ------------------------------------------------------------------------------

def morphological_filtering(sdata, filters):

  filename = sdata.path.stem  # 'blocco4_c26'

  # Use regex to extract blocco_key and samples_key
  match = re.match(r'(blocco\d+)_(\w+)', filename)
  if match:
      blocco_key, samples_key = match.group(1), match.group(2)
      
  else:
      raise ValueError(f"Could not parse blocco_key and samples_key from: {sdata_path}")

  # to add: way to extract blocco_key and then 
  try:
      # Assuming 'sdata' is your loaded SpatialData object
      features_df = features_extraction(sdata, nuclei_element_name = f"{blocco_key}_nuclei_boundaries")
      print("Successfully extracted features:")
  except ValueError as e:
      print(f"Error: {e}")
  
  # adding the morphological characteristics in the colData of our sdata object
  cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
  sdata['nuclei_counts_nop'].obs = sdata['nuclei_counts_nop'].obs.join(features_df[cols], how='left')

  filtered = nuclei_filtering(features_df, filters)

  # Get the cell_ids to keep
  cell_ids_to_keep = filtered.index.tolist()
  
  # Filter the GeoDataFrame by index
  sdata[f'{blocco_key}_filtered_nuclei'] = sdata[f'{blocco_key}_nuclei_boundaries'].loc[cell_ids_to_keep]
  
  return sdata

# versione copilot da verificare, la mia funziona comunque
# def morphological_filtering(sdata, filters):
#     """
#     Extracts blocco_key and sample_key from sdata path, computes morphological features,
#     filters nuclei based on user criteria, and updates sdata accordingly.
# 
#     Args:
#         sdata (SpatialData): Your SpatialData object.
#         filters (dict): Filtering criteria for nuclei features.
# 
#     Returns:
#         SpatialData: Updated sdata object.
#     """
# 
#     # --- 1. Extract keys from path ---
#     filename = sdata.path.stem  # e.g. 'blocco4_c26'
#     match = re.match(r'(blocco\d+)_(\w+)', filename)
#     if not match:
#         raise ValueError(f"Could not parse blocco_key and samples_key from: {sdata.path}")
#     blocco_key, samples_key = match.group(1), match.group(2)
# 
#     # --- 2. Feature extraction ---
#     try:
#         features_df = features_extraction(sdata, nuclei_element_name=f"{blocco_key}_nuclei_boundaries")
#         print("Successfully extracted features.")
#     except ValueError as e:
#         print(f"Error in features_extraction: {e}")
#         return sdata
# 
#     # --- 3. Filter nuclei features ---
#     filtered = nuclei_filtering(features_df, filters)
# 
#     # --- 4. Update filtered nuclei in sdata ---
#     cell_ids_to_keep = filtered.index.tolist()
#     # Make sure the index is compatible!
#     nuclei_gdf = sdata[f"{blocco_key}_nuclei_boundaries"]
#     sdata[f"{blocco_key}_filtered_nuclei"] = nuclei_gdf.loc[cell_ids_to_keep]
# 
#     # Optionally: update AnnData obs if needed (no redundant join)
#     # If you want to update .obs columns with new features:
#     # for col in ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']:
#     #     sdata['nuclei_counts_nop'].obs[col] = features_df[col]
# 
#     return sdata


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

# la mia brutta copia
# def features_extraction(sdata, nuclei_element_name = "nuclei_boundaries"):
#   '''
#   Function to extract some features from the nuclei shapes
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

# gemini 2.5 pro adjusted
def features_extraction(sdata, nuclei_element_name="nuclei_boundaries"):
    """
    Extracts morphological features from nuclei shapes in a SpatialData object.

    This function rasterizes nuclei polygons into a label mask, then computes
    shape properties (e.g., area, eccentricity) for each nucleus.

    Args:
        sdata (spatialdata.SpatialData): The SpatialData object containing the nuclei shapes.
        nuclei_element_name (str): The key for the nuclei shapes element in sdata.

    Returns:
        pd.DataFrame: A DataFrame containing morphological features for each nucleus,
                      indexed by the original cell ID.

    Raises:
        ValueError: If the required elements or attributes are not found in sdata.
    """
    # 1. --- Input Validation ---
    # if not sopa_attrs_check(sdata):
    #     raise ValueError("SOPA attributes check failed. Please verify the sdata object.")
    if nuclei_element_name not in sdata.shapes:
        raise ValueError(f"'{nuclei_element_name}' not found in sdata.shapes.")
    # 2. --- Extract Coordinate System Key ---
    try:
        # Assumes a pattern like 'blocco1' is the coordinate system name
        blocco_key_match = re.search(r'(blocco\d+)', sdata.attrs['boundaries_shapes'])
        if not blocco_key_match:
            raise ValueError("Could not find a 'bloccoN' key in sdata.attrs['boundaries_shapes'].")
        blocco_key = blocco_key_match.group(1)
    except (KeyError, AttributeError):
        raise ValueError("sdata.attrs['boundaries_shapes'] is missing or invalid.")
    # 3. --- Rasterize Nuclei Polygons ---
    raster_key = f"{blocco_key}_raster_nuclei"
    element_extent = sd.get_extent(sdata[nuclei_element_name], coordinate_system=blocco_key, exact=True)
    sdata[raster_key] = sd.rasterize(
        sdata[nuclei_element_name],
        axes=["x", "y"],
        min_coordinate=[element_extent['x'][0], element_extent['y'][0]],
        max_coordinate=[element_extent['x'][1], element_extent['y'][1]],
        target_coordinate_system=blocco_key,
        target_unit_to_pixels=1,
    )
    # 4. --- Prepare Label Mask ---
    # Squeeze the array to (H, W) and ensure it's an integer type
    label_mask = sdata[raster_key].values.squeeze().astype(np.int32)
    # 5. --- Compute Shape Features (regionprops) ---
    properties_to_extract = [
        'label', 'area', 'eccentricity', 'solidity', 'extent',
        'major_axis_length', 'minor_axis_length'
    ]
    props_df = pd.DataFrame(regionprops_table(label_mask, properties=properties_to_extract))
    # 6. --- Map Labels to Original Cell IDs ---
    try:
        label_to_id = sdata[raster_key].attrs['label_index_to_category']
        props_df['cell_id'] = props_df['label'].map(label_to_id)
        props_df = props_df.set_index('cell_id')
        props_df = props_df.drop(columns='label')
    except KeyError:
        print("Warning: Could not find 'label_index_to_category' mapping. Index will not be set to cell_id.")
    return props_df

# --- Example Usage ---
# try:
#     # Assuming 'sdata' is your loaded SpatialData object
#     features_df = features_extraction(sdata, nuclei_element_name="my_nuclei")
#     print("Successfully extracted features:")
#     print(features_df.head())
# except ValueError as e:
#     print(f"Error: {e}")

# ------------------------------------------------------------------------------

def find_extreme_observations(
    sdata: sd.SpatialData,
    table_key: str,
    feature_columns: List[str],
    n_extremes: int = 3
) -> pd.DataFrame:
    """
    Finds the top and bottom N extreme observations for a list of features in an AnnData table.
    This version correctly handles cases where the index name is also a column.

    Args:
        sdata: The SpatialData object containing the table.
        table_key: The key of the AnnData table in sdata (e.g., 'nuclei_counts').
        feature_columns: A list of column names in the .obs DataFrame to analyze.
        n_extremes: The number of top and bottom observations to find for each feature.

    Returns:
        A pandas DataFrame summarizing the extreme observations.
    """
    try:
        # Work on a copy to avoid modifying the original sdata
        obs_df = sdata[table_key].obs.copy()
    except KeyError:
        raise KeyError(f"Table with key '{table_key}' not found in the SpatialData object.")

    # --- FIX ---
    # If the index name also exists as a column, drop the column.
    # The index is the single source of truth for the ID.
    if obs_df.index.name is not None and obs_df.index.name in obs_df.columns:
        obs_df = obs_df.drop(columns=[obs_df.index.name])

    all_extremes_list = []

    for column in feature_columns:
        if column not in obs_df.columns:
            print(f"⚠️ Warning: Column '{column}' not found in '{table_key}.obs'. Skipping.")
            continue

        top_n = obs_df.nlargest(n_extremes, column)
        bottom_n = obs_df.nsmallest(n_extremes, column)
        
        top_n = top_n.assign(feature=column, extreme_type=f'Top {n_extremes}')
        bottom_n = bottom_n.assign(feature=column, extreme_type=f'Bottom {n_extremes}')
        
        all_extremes_list.extend([top_n, bottom_n])

    if not all_extremes_list:
        return pd.DataFrame()

    summary_df = pd.concat(all_extremes_list)

    # Now, reset_index will not create a duplicate column
    summary_df = summary_df.reset_index()

    # Reorder columns for better readability
    all_original_cols = [col for col in obs_df.columns if col not in feature_columns]
    
    # The index column is now called 'cell_id' (or whatever it was named)
    id_col_name = obs_df.index.name
    display_cols = ['feature', 'extreme_type', id_col_name] + feature_columns + all_original_cols
    
    final_cols = [col for col in display_cols if col in summary_df.columns]
    
    return summary_df[final_cols]


# --- How to use the corrected function ---
# cols_to_analyze = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']
# 
# # This will now produce a DataFrame with no duplicate columns
# extreme_nuclei_df = find_extreme_observations(
#     sdata=sdata,
#     table_key='nuclei_counts',
#     feature_columns=cols_to_analyze,
#     n_extremes=3
# )
# 
# # Verify the columns are clean
# print(extreme_nuclei_df.columns)

# ------------------------------------------------------------------------------

def sopa_attrs_check(sdata):
    """
    Validates that the SpatialData object has the required SOPA attributes.
    
    Checks if sdata.attrs contains the expected 'bins_table' and 'cell_segmentation_image'
    keys with appropriate values based on the sample/blocco ID.
    
    Parameters
    ----------
    sdata : SpatialData
        The SpatialData object to validate.
        
    Returns
    -------
    bool
        True if all required attributes are present and have correct values, False otherwise.
        
    Notes
    -----
    This function expects 'bins_table' to be 'filtered' and 'cell_segmentation_image'
    to follow the pattern '{blocco}_full_image'.
    """
    # Extract required attributes
    attrs = sdata.attrs
    
    # Check bins_table attribute
    if 'bins_table' not in attrs:
        print("Missing attribute: 'bins_table'")
        return False
    if attrs['bins_table'] != 'filtered':
        print(f"Attribute 'bins_table' has value '{attrs['bins_table']}', expected 'filtered'")
        return False
        
    # Check cell_segmentation_image attribute
    if 'cell_segmentation_image' not in attrs:
        print("Missing attribute: 'cell_segmentation_image'")
        return False
        
    # Extract blocco from the cell_segmentation_image value
    cell_seg_image = attrs['cell_segmentation_image']
    expected_suffix = "_full_image"
    
    # Verify the cell_segmentation_image follows the expected pattern
    if not cell_seg_image.endswith(expected_suffix):
        print(f"Attribute 'cell_segmentation_image' does not end with '{expected_suffix}'")
        return False
        
    # For additional validation, verify the blocco prefix is consistent across attributes
    blocco = cell_seg_image.replace(expected_suffix, "")
    
    # If 'boundaries_shapes' is present, check that it starts with the same blocco
    if 'boundaries_shapes' in attrs and not attrs['boundaries_shapes'].startswith(blocco):
        print(f"Warning: 'boundaries_shapes' value '{attrs['boundaries_shapes']}' doesn't match expected blocco '{blocco}'")
        # Not returning False here as you mentioned only checking the two main attributes
    
    # If 'tissue_segmentation_image' is present, check that it has the same value as cell_segmentation_image
    if 'tissue_segmentation_image' in attrs and attrs['tissue_segmentation_image'] != cell_seg_image:
        print(f"Warning: 'tissue_segmentation_image' value '{attrs['tissue_segmentation_image']}' doesn't match 'cell_segmentation_image'")
        # Not returning False here as you mentioned only checking the two main attributes
    
    return True

# ------------------------------------------------------------------------------

# Problem with the aggregation step, so let's try to divide the 2, segment all and 
# after that post process all.

# def apply_segmentation(sdata, expand_radius_ratio = None, no_overlap = True, filters = None, **kwargs):
  # '''
  # Function to apply segmentation to our data, one sample at a time
  # 
  # attrs structure:
  #   sdata.attrs:
  #     {'bins_table': 'filtered',
  #      'boundaries_shapes': 'blocco9_intissue',
  #      'cell_segmentation_image': 'blocco9_full_image',
  #      'tissue_segmentation_image': 'blocco9_full_image'}
  # '''
  # # Set-up keys to define all the elements needed
  # filename = sdata.path.stem  # 'blocco4_c26'
  # 
  # # Use regex to extract blocco_key and samples_key
  # match = re.match(r'(blocco\d+)_(\w+)', filename)
  # if match:
  #     blocco_key, samples_key = match.group(1), match.group(2)
  #     
  # else:
  #     raise ValueError(f"Could not parse blocco_key and samples_key from: {sdata_path}")
  # 
  # # set-up sopa metadata
  # sopa_attrs_check(sdata)
  # 
  # # 2. Segmentation
  # 
  # # # intissue polygons cleaned - added in the function to preprocess the samples
  # # sdata['blocco4_intissue'] = sdata['blocco4_intissue'][sdata[f'{blocco_key}_intissue'].name==f'{samples_key}']
  # 
  # # 2a. Divide the images in patches (overlapping if wanted)
  # sopa.make_image_patches(sdata, patch_width = 300, roi_key = f"{blocco_key}_intissue")
  # 
  # # 2c. Segmentation with stardist algorithm
  # sopa.segmentation.stardist(sdata, model_type='2D_versatile_he', 
  # image_key= sdata.attrs['cell_segmentation_image'], min_area=10, delete_cache=True, 
  # recover=False, prob_thresh=0.2, nms_thresh=0.6, key_added = f'{blocco_key}_nuclei_boundaries')
  # 
  # # 3. Post processing
  # 
  # # 3a. nuclei aggregation with no overlap of bins and nuclei after expansion
  # sopa.aggregate(sdata, key_added = 'nuclei_counts_nop', bins_key= "filtered", 
  # shapes_key = f"{blocco_key}_nuclei_boundaries", expand_radius_ratio=expand_radius_ratio, min_transcripts=1, 
  # min_intensity_ratio=0.1, no_overlap = no_overlap)
  # 
  # 
  # # 3b. Nuclei filtering based on morphological features
  # # if filters not defined, use this one
  # if filters is None:
  #   filters = {'area': (50, 4000),  # eliminiamo i più piccoli e quelli troppo grandi
  #   'eccentricity': (None, 0.95), # eliminiamo i nuclei più parabolici (di solito con un lato simil retto)
  #   'solidity': (0.7, None), # eliminiamo i più concavi
  #   'extent': (0.2, None)} # eliminiamo i nuclei più irregolari (un misto tra concavi e parabolici)
  #   sdata = morphological_filtering(sdata, filters)
  # else:
  #   sdata = morphological_filtering(sdata, filters)
  # 
  # # 3c. Annotating the table with the spatial element (nuclei polys)
  # sdata["nuclei_counts_nop"].obs["region"] = "blocco4_filtered_nuclei"
  # sdata.set_table_annotates_spatialelement("nuclei_counts_nop", region = f"{blocco_key}_filtered_nuclei", region_key="region", instance_key="cell_id")
  # 
  # # 3d. matching table with the filtered nuclei
  # sdata['nuclei_counts_nop'] = sd.match_table_to_element(sdata, element_name = f"{blocco_key}_filtered_nuclei", table_name='nuclei_counts_nop')
  # 
  # # 3e. Filtering bins_gdf 
  # sdata[f'{blocco_key}_intissue_002um']['cell_id'] = bin_to_cell_id_vector(sdata, table_key = 'nuclei_counts_nop')
  # sdata[f'{blocco_key}_intissue_filter'] = sdata[f'{blocco_key}_intissue_002um'][sdata[f'{blocco_key}_intissue_002um']['cell_id'].notna()]
  # 
  # # If you want to annotate the table of the genes vs nuclei with the filtered bins 
  # # sdata['nuclei_counts_nop'].obs['region'] = 'blocco4_intissue_filter'
  # # sdata.set_table_annotates_spatialelement('nuclei_counts_nop', region = 'blocco4_intissue_filter', region_key=None, instance_key = 'cell_id')
  # 
  # return sdata

def segmentation_step(sdata):
  '''
  Function to apply segmentation to our data, one sample at a time

  attrs structure:
    sdata.attrs:
      {'bins_table': 'filtered',
       'boundaries_shapes': 'blocco9_intissue',
       'cell_segmentation_image': 'blocco9_full_image',
       'tissue_segmentation_image': 'blocco9_full_image'}
  '''
  # Set-up keys to define all the elements needed
  filename = sdata.path.stem  # 'blocco4_c26'

  # Use regex to extract blocco_key and samples_key
  match = re.match(r'(blocco\d+)_(\w+)', filename)
  if match:
      blocco_key, samples_key = match.group(1), match.group(2)
      
  else:
      raise ValueError(f"Could not parse blocco_key and samples_key from: {sdata_path}")

  # set-up sopa metadata
  sopa_attrs_check(sdata)
  
  # 2. Segmentation
  
  # # intissue polygons cleaned - added in the function to preprocess the samples
  # sdata['blocco4_intissue'] = sdata['blocco4_intissue'][sdata['blocco4_intissue'].name=="c26"]
  
  # 2a. Divide the images in patches (overlapping if wanted)
  sopa.make_image_patches(sdata, patch_width = 300, roi_key = f"{blocco_key}_intissue")
  
  # 2c. Segmentation with stardist algorithm
  sopa.segmentation.stardist(sdata, model_type='2D_versatile_he', 
  image_key= sdata.attrs['cell_segmentation_image'], min_area=10, delete_cache=True, 
  recover=False, prob_thresh=0.2, nms_thresh=0.6, key_added = f'{blocco_key}_nuclei_boundaries')
  
  return sdata
  
# Post processing function 

def postprocess_step(sdata, expand_radius_ratio = None, no_overlap = True, filters = None, **kwargs):
  '''
  Function to apply post processing to our data, one sample at a time

  attrs structure:
    sdata.attrs:
      {'bins_table': 'filtered',
       'boundaries_shapes': 'blocco9_intissue',
       'cell_segmentation_image': 'blocco9_full_image',
       'tissue_segmentation_image': 'blocco9_full_image'}
  '''
  # Set-up keys to define all the elements needed
  filename = sdata.path.stem  # 'blocco4_c26'

  # Use regex to extract blocco_key and samples_key
  match = re.match(r'(blocco\d+)_(\w+)', filename)
  if match:
      blocco_key, samples_key = match.group(1), match.group(2)
      
  else:
      raise ValueError(f"Could not parse blocco_key and samples_key from: {sdata_path}")

  # set-up sopa metadata
  sopa_attrs_check(sdata)
  
  # 3. Post processing
  
  # 3a. nuclei aggregation with no overlap of bins and nuclei after expansion
  sopa.aggregate(sdata, key_added = 'nuclei_counts_nop', bins_key= "filtered", 
  shapes_key = f"{blocco_key}_nuclei_boundaries", expand_radius_ratio=expand_radius_ratio, min_transcripts=1, 
  min_intensity_ratio=0.1, no_overlap = no_overlap)
  
  
  # 3b. Nuclei filtering based on morphological features
  # if filters not defined, use this one
  if filters is None:
    filters = {'area': (50, 4000),  # filter out the bigger and smaller nuclei
    'eccentricity': (None, 0.95), # filter the most parabolic nuclei (one edge almost straight)
    'solidity': (0.7, None), # filter the most concave
    'extent': (0.2, None)} # filter the most irregular nuclei (mix between parabolic and concave)
    sdata = morphological_filtering(sdata, filters)
  else:
    sdata = morphological_filtering(sdata, filters)
  
  # 3c. Annotating the table with the spatial element (nuclei polys)
  sdata["nuclei_counts_nop"].obs["region"] = f"{blocco_key}_filtered_nuclei"
  sdata.set_table_annotates_spatialelement("nuclei_counts_nop", region = f"{blocco_key}_filtered_nuclei", region_key="region", instance_key="cell_id")
  
  # 3d. matching table with the filtered nuclei
  sdata['nuclei_counts_nop'] = sd.match_table_to_element(sdata, element_name = f"{blocco_key}_filtered_nuclei", table_name='nuclei_counts_nop')
  
  # 3e. Filtering bins_gdf 
  sdata[f'{blocco_key}_intissue_002um']['cell_id'] = bin_to_cell_id_vector(sdata, table_key = 'nuclei_counts_nop')
  sdata[f'{blocco_key}_intissue_filter'] = sdata[f'{blocco_key}_intissue_002um'][sdata[f'{blocco_key}_intissue_002um']['cell_id'].notna()]
  
  # If you want to annotate the table of the genes vs nuclei with the filtered bins 
  # sdata['nuclei_counts_nop'].obs['region'] = 'blocco4_intissue_filter'
  # sdata.set_table_annotates_spatialelement('nuclei_counts_nop', region = 'blocco4_intissue_filter', region_key=None, instance_key = 'cell_id')

  return sdata

































