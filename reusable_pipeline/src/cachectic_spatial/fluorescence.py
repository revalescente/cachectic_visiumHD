import sopa


def add_fluorescence_to_table(
    sdata,
    image_key,
    shapes_key,
    table_key,
    column_name="GFP_value",
):
    values = sopa.aggregation.aggregate_channels(
        sdata,
        image_key=image_key,
        shapes_key=shapes_key,
        expand_radius_ratio=0,
        mode="max",
        no_overlap=False,
    )

    sdata.tables[table_key].obs[column_name] = values.max(axis=1)
