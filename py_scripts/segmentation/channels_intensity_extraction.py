import spatialdata as sd
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Identity, get_transformation, remove_transformation, Translation, set_transformation
import sopa
import spatialdata_plot
import matplotlib.pyplot as plt
from skimage import io
import numpy as np
from py_scripts.utils.utils_fun import read_from_json
from spatialdata_plot.pl.utils import set_zero_in_cmap_to_transparent

samples_dict = read_from_json("/mnt/europa/valerio/repositories/cachetic_visiumHD/json/blocco_sample_bbox_dict.json")
samples_key = [
    details['sample_key'] 
    for blocco in samples_dict.values() 
    for details in blocco.values()
]

sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/samples/blocco1_sham.zarr")

img = sdata['full_image'].values()

def staining(self) -> np.ndarray:
    if self.channel is None:
        thumbnail_2d = np.array(self.image.max(dim="c").transpose("y", "x"))
    else:
        thumbnail_2d = np.array(self.image.sel(c=self.channel).transpose("y", "x")
    threshold = np.quantile(thumbnail_2d, self.clip_parameters[0]) / self.clip_parameters[1]
    thumbnail_2d = ((thumbnail_2d / threshold).clip(0, 1) * 255).astype(np.uint8)
