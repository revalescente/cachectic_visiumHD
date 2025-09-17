sd.get_extent(sdata['blocco4_nuclei_boundaries'], coordinate_system='blocco4', exact=True)
# {'x': (1627.5128723751084, 11474.5), 'y': (9293.514379291657, 15990.5)}


# rasterize shapes
sdata["raster_nuclei"] = sd.rasterize(
    sdata["blocco4_nuclei_boundaries"],
    ["x", "y"],
    min_coordinate=[1627.5128723751084, 9293.514379291657],
    max_coordinate=[11474.5, 15990.5],
    target_coordinate_system="blocco4",
    target_unit_to_pixels=1,
)

# transform into integer values

# If your array is an xarray with shape (1, H, W), get the first channel
label_mask = sdata['raster_nuclei'].values

# If shape is (1, H, W), squeeze to (H, W)
if label_mask.ndim == 3 and label_mask.shape[0] == 1:
    label_mask = label_mask[0]

# Cast to integer type
label_mask = label_mask.astype(np.int32)

# add features of the rasterized nuclei boundaries

# 2. Compute regionprops
props = regionprops_table(label_mask, properties=[
    'label', 'area', 'eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length'])
props_df = pd.DataFrame(props)

# Get mapping dictionary from xarray attributes
label_to_id = sdata['raster_nuclei'].attrs['label_index_to_category']

# Map: create a new column with the nucleus ID
props_df['cell_id'] = props_df['label'].map(label_to_id)
props_df = props_df.drop(columns='label')
# props_df.to_csv('/mnt/europa/valerio/nuclei_features_df.csv', index=False)

# DA VERIFICARE CHE LA FUNZIONE SCRITTA DA GEMINI SIA BUONA! E FUNZIONI
# CONFRONTARE I DUE RISULTATI SE SONO UGUALI

try:
    # Assuming 'sdata' is your loaded SpatialData object
    features_df = features_extraction(sdata, nuclei_element_name="blocco4_nuclei_boundaries")
    print("Successfully extracted features:")
    print(features_df.head())
except ValueError as e:
    print(f"Error: {e}")

# testing equality
props_df.equals(features_df)
# correct


# ------------------------------------------------------------------------------

# I'd like to make a couple plot of the most specific nuclei
# i need to find the nuclei with the specific features values 
# get the extent, bounding box and then plot them? i guess so

# extract the cell_id of the nuclei to plot:
cols = ['eccentricity', 'solidity', 'extent', 'major_axis_length', 'minor_axis_length']

# A list to hold all the results
all_extremes = []

# Loop through each feature and get the extremes
for column in cols:
    if column in sdata['nuclei_counts'].obs.columns:
        top_3, bottom_3 = get_extreme_observations(sdata['nuclei_counts'].obs, column, number_observations=3)
        
        # Add a column to identify the feature and extreme type
        top_3 = top_3.assign(feature=column, extreme_type='Top 3')
        bottom_3 = bottom_3.assign(feature=column, extreme_type='Bottom 3')
        
        all_extremes.extend([top_3, bottom_3])

# Combine everything into a single DataFrame for a nice summary
summary_df = pd.concat(all_extremes)

# Reorder columns for better readability
display_cols = ['feature', 'extreme_type', 'cell_id'] + feature_columns
summary_df = summary_df.reset_index().rename(columns={'index': 'cell_id'})
summary_df = summary_df[display_cols]


# Display the final summary table
print(summary_df.to_string())

# we have the info about the cell_id of the extremes nuclei. 
# now let's try to plot one nucleus

# 1. Define the list of cell_ids you want to keep.
#    Replace these example IDs with your actual list.
ids_to_keep = summary_df.index

# 2. Filter the GeoDataFrame using .loc
#    This selects only the rows whose index (cell_id) is in your list.
nucleo_geo = sdata['blocco4_nuclei_boundaries'].loc[[summary_df.index[0]]]
sdata['nucleo_shape'] = nucleo_geo

sd.get_extent(sdata['nucleo_shape'], coordinate_system='blocco4', exact=True)
# {'x': (3955.5047917988886, 3966.4837517193064), 'y': (12629.535219786241, 12649.252859600709)}

plt.figure(figsize=(20, 20))
ax = plt.gca()
sdata.query.bounding_box(
    axes=["x", "y"],
    min_coordinate=[3955.5047917988886-50, 12629.535219786241-50],
    max_coordinate=[3966.4837517193064+50, 12649.252859600709+50],
    target_coordinate_system="blocco4",
).pl.render_images("blocco4_full_image").pl.render_shapes("blocco4_nuclei_boundaries", outline=True, outline_alpha=1, outline_width=1.5, fill_alpha=0
).pl.render_shapes("nucleo_shape", color = "red", fill_alpha = 0.5).pl.show(ax = ax, coordinate_systems="blocco4", save = 'output_python/nuclei_Strange/most_ecc_1.png')










