from dataclasses import dataclass


@dataclass(frozen=True)
class Paths:
    samples_json: str
    visium: str
    fullres_image: str
    annotations: str
    fluorescence_image: str
    output: str


@dataclass(frozen=True)
class SegmentedObject:
    name: str
    segmentation_path: str
    labels_key: str
    shapes_key: str
    table_key: str
    raster_key: str
    morphology_properties: tuple[str, ...]
    min_transcripts: int = 10


PATHS = Paths(
    samples_json="../json/blocco_sample_bbox_dict.json",
    visium="/mnt/europa/valerio/data/visium/{block}/outs",
    fullres_image="/mnt/europa/valerio/HE_images/color_corrected/blocchi/pp_{block}_20x.tif",
    annotations="/mnt/europa/valerio/data/json/geojson_dir/intissue_GFP_polys/{block}_{sample}.geojson",
    fluorescence_image="/mnt/europa/valerio/Fluo_images/warped_tif/samples/{block}_{sample}_fluo_image.tif",
    output="/mnt/europa/valerio/data/zarr_store/objects/{block}_{sample}",
)

MORPHOLOGY_PROPERTIES = (
    "area",
    "eccentricity",
    "solidity",
    "extent",
    "major_axis_length",
    "minor_axis_length",
)

NUCLEI = SegmentedObject(
    name="nuclei",
    segmentation_path="/mnt/europa/valerio/data/arivis_cloud_segmentation/segmentation_masks/{block}_{sample}_finalprediction.tiff",
    labels_key="nuclei_labels",
    shapes_key="nuclei_shapes",
    table_key="nuclei_table",
    raster_key="nuclei_raster",
    morphology_properties=MORPHOLOGY_PROPERTIES,
)

FIBRES = SegmentedObject(
    name="fibres",
    segmentation_path="/mnt/europa/valerio/data/fibre_segmentation/{block}_{sample}_fibres.tiff",
    labels_key="fibres_labels",
    shapes_key="fibres_shapes",
    table_key="fibres_table",
    raster_key="fibres_raster",
    morphology_properties=MORPHOLOGY_PROPERTIES,
)

SEGMENTED_OBJECTS = (NUCLEI, FIBRES)
BIN_SIZES = ("016", "008", "002")
BLOCK_TO_PROCESS = "blocco1"
