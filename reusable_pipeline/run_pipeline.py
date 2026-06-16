from cachectic_spatial.block import process_block
from cachectic_spatial.config import BLOCK_TO_PROCESS, PATHS
from cachectic_spatial.input import read_samples


samples = read_samples(PATHS.samples_json)
process_block(BLOCK_TO_PROCESS, samples[BLOCK_TO_PROCESS])
