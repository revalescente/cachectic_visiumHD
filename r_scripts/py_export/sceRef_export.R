library(zellkonverter)
library(SingleCellExperiment)

sce_ref <- readRDS(ref_path)

writeH5AD(sce_ref, file = paste0(ref_path, "multiome_anndata.h5ad"), X_name = "counts", verbose = TRUE)