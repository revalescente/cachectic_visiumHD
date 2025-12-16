library(SpatialExperiment)
library(arrow)
library(ggplot2)
library(dplyr)

load("~/data/Rdata/spe_blocco1_emma.RData")
blocco1

blocco1$c26STAT3_b1

parquet_path = "/mnt/europa/data/sandri/241219_A00626_0902_AHWH77DMXY_3/space_out_3.1/blocco1/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet"
bin_poly <- read_parquet(parquet_path)

head(bin_poly)
head(colnames(blocco1$sham_b1))

summary(bin_poly$pxl_row_in_fullres)
summary(bin_poly$pxl_col_in_fullres)
