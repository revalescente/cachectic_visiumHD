# BiocManager::install("HelenaLC/SpatialData")
# BiocManager::install("keller-mark/anndataR", ref="spatialdata")
# BiocManager::install("keller-mark/pizzarr")

# let's try with Rarr
# install.packages("Rarr", repos = "https://bioconductor.org/packages/release/bioc")

# some library need to be updated, let's try to install locally
install.packages(c("stringfish", "RApiSerialize", "qs"), lib="R/x86_64-pc-linux-gnu-library/4.5")
library(ggplot2)
library(ggnewscale)
library(SpatialData)
library(SingleCellExperiment)

readSpatialData


(b1_stat3 <- readSpatialData("data/b1_stat3.zarr", anndataR=TRUE))
