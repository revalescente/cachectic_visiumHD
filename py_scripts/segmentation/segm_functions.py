import geopandas as gpd
import pandas as pd
from spatialdata.models import ShapesModel

# Function to spatial join bins with respect to the nuclei
def assign_bins_to_nearest_nucleus(bins, nuclei_bounds):
    """
    For each bin that intersects one or more nuclei, assign it to the nearest nucleus boundary
    (by centroid-to-boundary distance), and return a GeoDataFrame with bin-to-nucleus assignment.
    """
    # Ensure location_id is categorical for efficient groupby (optional)
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
# filters = {'area': (1000, None), 'eccentricity': (None, 0.8), 'solidity': (0.95, None)}
# filtered = nuclei_filtering(props_df, filters)
