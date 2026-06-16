import geopandas as gpd
import spatialdata as sd


def annotate_table(
    sdata,
    table_key,
    shapes_key,
    annotations_key,
    coordinate_system,
    sample,
    filter_in_tissue=False,
):
    shapes = sdata.transform_element_to_coordinate_system(shapes_key, coordinate_system)
    annotations = sdata.transform_element_to_coordinate_system(
        annotations_key,
        coordinate_system,
    )
    joined = gpd.sjoin(shapes, annotations, how="inner", predicate="intersects")

    table = sdata.tables[table_key]
    instance_key = table.uns["spatialdata_attrs"]["instance_key"]

    table.obs["in_treatment"] = False
    table.obs["to_discard"] = False
    table.obs["in_tissue"] = False

    bins_by_annotation = joined.groupby("name").apply(lambda rows: rows.index.tolist())

    for annotation_name, bin_ids in bins_by_annotation.items():
        column = annotation_column(annotation_name, sample)
        if column:
            table.obs.loc[table.obs[instance_key].isin(bin_ids), column] = True

    sdata.tables[table_key] = table

    if filter_in_tissue:
        sdata.tables[table_key] = table[table.obs["in_tissue"]].copy()
        matched = sd.match_element_to_table(sdata, shapes_key, table_key)
        sdata.shapes[shapes_key] = matched[0][shapes_key]


def annotation_column(annotation_name, sample):
    if "fibre_trattate" in annotation_name:
        return "in_treatment"
    if "infiammazione" in annotation_name:
        return "to_discard"
    if sample in annotation_name:
        return "in_tissue"
