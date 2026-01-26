import geopandas as gpd
import spatialdata as sd

# nuclei importati da qupath
geojson_path = "/mnt/europa/valerio/data/json/geojson_dir/annotated_nuclei/nuclei_b1_sham.geojson"
nuclei_poly = gpd.read_file(geojson_path)
nuclei_poly = nuclei_poly.set_crs(None, allow_override=True)

# non vanno bene perche perdono l'informazione dell' "cell_id", devo fare in modo che si salvi:
# 1- provo a inserire una colonna e vedere se rimane, dubito
# 2- provo a inserire il cell_id come "id" univoco del oggetto
# 3- provo a inserire nella colonna "name"
# 4- idee esaurite

spatialdata_path =  "/mnt/europa/valerio/data/zarr_store/samples/blocco1_sham.zarr"
sdata = sd.read_zarr(spatialdata_path)

shapes_dict, _ = sd.match_element_to_table(sdata, "nuclei_boundaries", "nuclei_counts")
gpd = shapes_dict['nuclei_boundaries']

gpd['cell_id'] = gpd.index
gpd['name'] = gpd.index

# gdf = gdf.reset_index()
# 5. Save to disk
# Parquet is the best format for GeoDataFrames (fast, preserves types)
output_dir = "/mnt/europa/valerio/data/extracted_shapes"
output_path = os.path.join(output_dir, "blocco1_sham.geojson")

gpd.to_file(output_path, driver="GeoJSON")
