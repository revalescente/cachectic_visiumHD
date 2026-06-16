from cachectic_spatial.config import PATHS
from cachectic_spatial.input import read_visium
from cachectic_spatial.sample import process_sample


def process_block(block, samples):
    sdata = read_visium(block, PATHS)

    for sample, sample_info in samples.items():
        process_sample(sdata, block, sample, sample_info)
