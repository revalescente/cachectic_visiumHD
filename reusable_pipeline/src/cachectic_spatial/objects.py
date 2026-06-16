import sopa

from cachectic_spatial.annotations import annotate_table
from cachectic_spatial.elements import add_segmented_shapes
from cachectic_spatial.fluorescence import add_fluorescence_to_table
from cachectic_spatial.input import read_segmentation
from cachectic_spatial.morphology import add_morphology_to_table, extract_morphology
from cachectic_spatial.transcripts import aggregate_transcripts


def add_segmented_objects(
    sdata,
    block,
    sample,
    object_spec,
    transformation,
    fluorescence_key,
    annotations_key,
    bins_key="square_002um",
):
    mask = read_segmentation(block, sample, object_spec)

    add_segmented_shapes(
        sdata,
        mask,
        object_spec.labels_key,
        object_spec.shapes_key,
        block,
        transformation,
    )

    sopa.utils.set_sopa_attrs(
        sdata,
        cell_segmentation_key=f"{block}_full_image",
        tissue_segmentation_key=f"{block}_full_image",
        transcripts_key=None,
        boundaries_key=object_spec.shapes_key,
        bins_table_key=bins_key,
    )

    aggregate_transcripts(
        sdata,
        object_spec.shapes_key,
        bins_key,
        object_spec.table_key,
        object_spec.min_transcripts,
    )

    add_fluorescence_to_table(
        sdata,
        fluorescence_key,
        object_spec.shapes_key,
        object_spec.table_key,
    )

    morphology = extract_morphology(
        sdata,
        object_spec.shapes_key,
        object_spec.raster_key,
        block,
        object_spec.morphology_properties,
    )
    add_morphology_to_table(sdata, object_spec.table_key, morphology)

    table = sdata.tables[object_spec.table_key]
    table.obs["object_type"] = object_spec.name
    table.obs["block_id"] = block
    table.obs["sample_id"] = sample

    annotate_table(
        sdata,
        object_spec.table_key,
        object_spec.shapes_key,
        annotations_key,
        block,
        sample,
    )
