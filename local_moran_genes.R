library(SpatialExperiment)
library(arrow)
library(sf)
library(dplyr) # Optional, for data manipulation
library(SpatialFeatureExperiment)
library(Voyager)

# spatialexperiment of emma binned data - blocco 1 sham of interest
spe_emma <- readRDS("/mnt/europa/valerio/data/Rdata/emma_16um_bin_annotated_spelist.RDS")

spe <- spe_emma$c26_b4

# read the parquet of the bins
bins <- read_parquet("/mnt/europa/valerio/data/zarr_store/binned_samples/version_1.0.0/blocco4_c26/shapes/blocco4_square_016um/shapes.parquet")
bins$geometry <- st_as_sfc(bins$geometry, crs = NA) 
bins_sf <- st_as_sf(bins)
print(bins_sf)


# create the SpatialFeatureExperiment
spe

colData(spe)

plot(st_geometry(bins_sf))

# read zarr file because for some reason shapes parquet doesn't have the right index 
# sdata = sd.read_zarr("/mnt/europa/valerio/data/zarr_store/binned_samples/version_1.0.0/blocco1_sham")

# spatial join ----
spe_coords <- spatialCoords(spe) |> as.data.frame()
spe_points_sf <- st_as_sf(spe_coords, coords = c("pxl_col_in_fullres", "pxl_row_in_fullres"), crs = st_crs(bins_sf))
spe_points_sf$bin_id <- colnames(spe)

# st_join defaults to "intersects" (point inside polygon)
joined_data <- st_join(bins_sf, spe_points_sf, join = st_intersects, left = FALSE)

# Check the result
head(joined_data)

matched_info <- joined_data[match(colnames(spe), joined_data$bin_id), ]
head(matched_info)

# creation of the SpatialFeatureExperiment object
imgData(spe) <- NULL
sfe <- toSpatialFeatureExperiment(spe)

# sfe <- SpatialFeatureExperiment(assays = assays(spe),
#                                 colData = colData(spe),
#                                 #spatialCoordsNames = c("pxl_col_in_fullres", "pxl_row_in_fullres"),
#                                 #spatialCoords = spatialCoords(spe),
#                                 #colGeometries = list("bins_polygons" = matched_info),
#                                 #spotDiameter = 2,
#                                 #unit = "full_res_image_pixel"
#                                 )

sfe

# lognormcounts
sfe <- scuttle::logNormCounts(sfe)

# features of interest, for example 
features <- c("Trim63", "Myh4")
colGraphName <- "knn8"
colGeometryName <- "bins_polygons"
pointsize <- 2

# not working and even useless
Voyager::plotSpatialFeature(sfe, features,
                   colGeometryName = "bins_polygons", 
                   nrow = 2, size = pointsize, scattermore = FALSE) +
  theme_void()

# neighborhood graphs for the bins
colGraph(sfe, "knn8") <-
  findSpatialNeighbors(sfe,
                       method = "knearneigh", # wraps the spdep function with the same name
                       k = 8,
                       zero.policy = TRUE
  )

Voyager::plotColGraph(sfe,
                   colGraphName = "knn8",
                   colGeometryName = "centroids",
                   bbox = c(xmin = 6000, xmax = 8000, ymin = 10000, ymax = 11000)
                   
)

# local Moran's I ----
sfe <- Voyager::runUnivariate(sfe,
                     features = features[1],
                     colGraphName = "knn8",
                     type = "localmoran")

sfe <- Voyager::runBivariate(sfe,
                              feature1 = features[1],
                              feature2 = features[2],
                              colGraphName = "knn8",
                              type = "localmoran_bv",
                              nsim = 100)

p_f2_f1 <- Voyager::plotLocalResult(sfe, 
                         name = "localmoran_bv",       # The method you used
                         features = paste(features[2], features[1], sep="__"),    # The specific pair you calculated
                         colGeometryName = "centroids", # Your polygon geometry name
                         divergent = TRUE, 
                         diverge_center = 0,
                         size = 1)

p_f1_f2 <- Voyager::plotLocalResult(sfe, 
                                    name = "localmoran_bv",       # The method you used
                                    features = paste(features[1], features[2], sep="__"),    # The specific pair you calculated
                                    colGeometryName = "centroids", # Your polygon geometry name
                                    divergent = TRUE, 
                                    diverge_center = 0,
                                    size = 1)

# almost simmetrical but not really

sfe <- runUnivariate(
  sfe,
  features[1],
  colGraphName = "knn8",
  type = "moran.plot",
  nsim = 400
)
sfe <- runUnivariate(
  sfe,
  features[2],
  colGraphName = "knn8",
  type = "moran.plot",
  nsim = 400
)

moranPlot(sfe,
          features[2],
          graphName = "knn8",
          swap_rownames = "symbol")
moranPlot(sfe,
          features[1],
          graphName = "knn8",
          swap_rownames = "symbol")


sfe <- runUnivariate(
  sfe,
  features[1],
  colGraphName = "knn8",
  type = "localmoran"
)
sfe <- runUnivariate(
  sfe,
  features[2],
  colGraphName = "knn8",
  type = "localmoran"
)

p_l2 <- plotLocalResult(
  sfe,
  name = "localmoran",
  features = features[2],
  attribute = "mean",
  colGeometryName = "centroids",
  size = 1
)

p_l1 <- plotLocalResult(
  sfe,
  name = "localmoran",
  features = features[1],
  attribute = "mean",
  colGeometryName = "centroids",
  size = 1
)

localResults(sfe)[["localmoran"]][[1]] <- 
  localResults(sfe)[["localmoran"]][[1]] |> 
  mutate(locClust = ifelse(`-log10p_adj` > -log10(0.05), as.character(mean), "non-siginificant"))

localResults(sfe)[["localmoran"]][[2]] <- 
  localResults(sfe)[["localmoran"]][[2]] |> 
  mutate(locClust = ifelse(`-log10p_adj` > -log10(0.05), as.character(mean), "non-siginificant"))

p_l1 <- plotLocalResult(
  sfe,
  name = "localmoran",
  features = features[1],
  attribute = "locClust",
  colGeometryName = "centroids",
  size = 1
)
p_l2 <- plotLocalResult(
  sfe,
  name = "localmoran",
  features = features[2],
  attribute = "locClust",
  colGeometryName = "centroids",
  size = 1
)

