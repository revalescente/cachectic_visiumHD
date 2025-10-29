import re
import json
import pandas as pd
import geopandas as gpd
import sopa
import spatialdata as sd
from spatialdata.transformations import Identity
from spatialdata.models import (ShapesModel, TableModel, Image2DModel)
from py_scripts.utils.utils_fun import read_from_json

# Fixed cell_id format - as it is in maps_bins
def format_cell_id(cell_id):
    # Convert to string, pad with leading zeros to 9 digits, and add prefix/suffix
    return f"cellid_{cell_id:09d}-1"

# Import and prepare the single sample importation in sdata 
def geojson_preparation(geojson_path, samples_dict, sample_key, blocco_key):
  # read geojson
  sranger_nuclei = gpd.read_file(geojson_path)
  sranger_nuclei = sranger_nuclei.set_crs(None, allow_override=True)
  sranger_nuclei['area'] = sranger_nuclei.geometry.area
  # sdata preparation
  sranger_parsed = ShapesModel.parse(sranger_nuclei, transformations = {'blocco1': Identity()})
  temp_sd = sd.SpatialData(shapes={"sranger_nuclei": sranger_parsed})
  ## cutting
  min_coord = samples_dict[blocco_key][sample_key]['min_coordinate']
  max_coord = samples_dict[blocco_key][sample_key]['max_coordinate']
  sdata = temp_sd.query.bounding_box(
                axes=['x', 'y'],
                min_coordinate=min_coord,
                max_coordinate=max_coord,
                target_coordinate_system=blocco_key
  )
  # fix cell_id naming - lazy rangers
  sdata['sranger_nuclei'] = sdata['sranger_nuclei'].rename(columns={'cell_id': 'orig_cell_id'})
  sdata['sranger_nuclei']['cell_id'] = sdata['sranger_nuclei']['orig_cell_id'].apply(format_cell_id)
  return(sdata['sranger_nuclei'])


# Import and prepare the anndata with the maps_bins information
def table_integration(sranger_nuclei_gpd, nuclei_matrix_path, region_key = 'sranger_nuclei'):
  #
  adata = sc.read_10x_mtx(
      nuclei_matrix_path,  # Directory containing the files
      var_names='gene_symbols',  # Use gene symbols as variable names
      cache=True  # Cache the result for faster reloading
  )
  adata.obs['region'] = region_key
  adata.obs['cell_id'] = adata.obs.index
  #filter by nuclei
  #
  cell_ids_to_keep = set(sranger_nuclei['cell_id'].unique())
  print(f"Number of unique cell IDs in sranger_nuclei: {len(cell_ids_to_keep)}")
  #
  # Create a boolean mask for observations in adata that match the cell IDs
  mask = adata.obs.index.isin(cell_ids_to_keep)
  #
  # Filter the adata object
  adata_filtered = adata[mask].copy()
  return(TableModel.parse(adata_filtered))

