from skimage.color import rgb2gray
from skimage.filters import threshold_otsu
from skimage.morphology import closing, disk, remove_small_holes, remove_small_objects


def background_mask(image):
    gray = rgb2gray(image)
    threshold = threshold_otsu(gray)

    tissue = gray < threshold
    tissue &= image.max(axis=2) > 20
    tissue = closing(tissue, disk(3))
    tissue = remove_small_objects(tissue, max_size=100)
    tissue = closing(tissue, disk(3))
    tissue = remove_small_holes(tissue, max_size=10000)

    return ~tissue
