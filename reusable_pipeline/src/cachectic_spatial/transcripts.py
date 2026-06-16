import sopa


def aggregate_transcripts(
    sdata,
    shapes_key,
    bins_key,
    table_key,
    min_transcripts,
):
    sopa.aggregate(
        sdata,
        key_added=table_key,
        bins_key=bins_key,
        shapes_key=shapes_key,
        expand_radius_ratio=0,
        min_transcripts=min_transcripts,
        min_intensity_ratio=0.15,
        no_overlap=True,
    )
