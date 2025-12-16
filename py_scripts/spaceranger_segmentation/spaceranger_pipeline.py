# required packages
# pip install anndata==0.12.0 pydeseq2==0.5.2 squidpy==1.6.5 spatialdata[extra]==0.4.0 geosketch==1.3 harmonypy==0.0.10 igraph==0.11.8

import spatialdata as spd
import gc
import py_scripts.utils.spaceranger_utility as sr
from py_scripts.utils.utils_fun import read_from_json

"""# **Section 3: Conversion of Space Ranger Output to Zarr Format and SpatialData Object Creation**

We begin this section by converting the standard Visium HD output into a Zarr file, 
as the `spatialdata` library expects data in this format. If you are running this code locally, 
this code needs to be run only once, as the datasets can be loaded directly from the saved Zarr files afterward.

In the following code snippet, a dictionary is created where each sample name serves as a unique key. 
The value associated with each key is a list containing:
* the path and filename of the filtered feature-cell matrix in `h5` format.
* the location and name of the image to be stored in the `SpatialData` object.
* the scale factors JSON file, so that the cell segmentation results can beoverlaid onto the tissue image.
* the cell segmentation GeoJSON file, so the cell segmentation results can be visualized on the tissue image.
* the desired name for the Zarr file.

Each key-value pair in this dictionary is then processed by the `create_zarr` helper function. 
For this specific example, we only use the `tissue_hires_image.png`. 
Other images, such as a high-resolution microscope image or CytAssist image, can be added to the `SpatialData` object.

"""
# /mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs
# /mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/Images/HE

# blocco1 and blocco2 no cell expansion
# /mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0nocell/ 

# Create and save Zarr files for the cell segmentation outputs.
# sample = {
#   "blocco1":["/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/raw_feature_cell_matrix.h5",
#              "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/spatial/tissue_hires_image.png",
#              "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/spatial/scalefactors_json.json",
#              "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0/blocco1/outs/segmented_outputs/nucleus_segmentations.geojson",
#              "/mnt/europa/valerio/data/json/geojson_dir/tissue_hires_image_blocco1.geojson",
#              "blocco1_cell"]
# }
# no cell
sample = {
  "blocco1":["/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco1/outs/segmented_outputs/raw_feature_cell_matrix.h5",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco1/outs/segmented_outputs/spatial/tissue_hires_image.png",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco1/outs/segmented_outputs/spatial/scalefactors_json.json",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco1/outs/segmented_outputs/nucleus_segmentations.geojson",
             "/mnt/europa/valerio/data/json/geojson_dir/tissue_hires_image_blocco1.geojson",
             "blocco1",
             "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/cell_expans_2um/blocchi/"
             ],
  "blocco2":["/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco2/outs/segmented_outputs/raw_feature_cell_matrix.h5",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco2/outs/segmented_outputs/spatial/tissue_hires_image.png",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco2/outs/segmented_outputs/spatial/scalefactors_json.json",
             "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_4.0_2micron/blocco2/outs/segmented_outputs/nucleus_segmentations.geojson",
             "/mnt/europa/valerio/data/json/geojson_dir/tissue_hires_image_blocco2.geojson",
             "blocco2",
             "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/cell_expans_2um/blocchi/"
            ]
}

print("Saving zarr files")
for key, inputs in sample.items():
    sr.create_zarr(count_matrix_path=inputs[0],
                image_path=inputs[1],
                scale_factors_path=inputs[2],
                geojson_path=inputs[3],
                intissue_geojson_path=inputs[4],
                sample_name=inputs[5],
                save_path = inputs[6])

del sample, inputs, key
gc.collect()

"""
Let's divide the samples in the blocco, I need a samples dictionary to know where to cut the bbox
"""
save_path = "/mnt/europa/valerio/data/zarr_store/spaceranger_v4/try2/"

blocchi = ["blocco1","blocco2"]

for blocco in blocchi:
  sdata = spd.read_zarr(f"/mnt/europa/valerio/data/zarr_store/spaceranger_v4/no_cell_expans/blocchi/{blocco}")
  sr.query_and_crop_better(sdata, save_path)

