# BiocManager::install("HelenaLC/SpatialData")
# BiocManager::install("keller-mark/anndataR", ref="spatialdata")
# BiocManager::install("keller-mark/pizzarr")
# BiocManager::install("HelenaLC/SpatialData.plot")

# let's try with Rarr
# install.packages("Rarr", repos = "https://bioconductor.org/packages/release/bioc")

# some library need to be updated, let's try to install locally
# install.packages(c("stringfish", "RApiSerialize", "qs"), lib="R/x86_64-pc-linux-gnu-library/4.5")
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(SpatialData)
library(SpatialData.plot)
library(SingleCellExperiment)


# Read data 
#(sdata2 <- readSpatialData("/mnt/europa/valerio/data/zarr_store/blocchi/blocco3_sham.zarr", anndataR=TRUE))
#(sdata <- readSpatialData("/mnt/europa/valerio/data/zarr_store/general/b1_stat3.zarr", anndataR=TRUE))
(sdata <- readSpatialData("/mnt/europa/valerio/data/zarr_store/concat4samples.zarr", anndataR=TRUE))

shapes(sdata2)
data(shape(sdata))$metadata
metadata(shape(sdata2, 'blocco4_nuclei_boundaries'))

# tables
hasTable(sdata, "blocco2_c26_filtered_nuclei")
# retrieve 'table' for an element
getTable(sdata, i <- "blocco2_c26murf1_filtered_nuclei")

(spe <- tables(sdata)$nuclei_counts_nop)
class(spe)
meta(tables(sdata2)$nuclei_counts_nop)

fim <- image(sdata2, "blocco4_full_image")

# different scales with different resolution
vapply(fim@data, dim, numeric(3))




