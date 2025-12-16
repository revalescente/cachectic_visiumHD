library(scuttle)
library(scater)
library(edgeR)
library(scran)
library(SpatialExperiment)
library(patchwork)
library(ggplot2)
library(SpaceTrooper)

Rd_path = "/mnt/europa/valerio/data/Rdata/"
# pb_ct <- readRDS(paste0(Rd_path, "pseudobulk_celltype_SingleR"))
spe_list <- readRDS(paste0(Rd_path, "spe_list_nuclei_ann"))

# pseudo bulk ----

dec.pbmc <- modelGeneVar(spe_prova)

# Visualizing the fit:
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)


# sce o spe ----

spe_prova <- spe_list$blocco2_c26
plotExpression(spe_prova, x="labels_singler", features="Ttn", colour="labels_singler")


spe_prova$library_size <- librarySizeFactors(spe_prova)
spe_prova <- logNormCounts(spe_prova, size.factors=spe_prova$library_size)
assayNames(spe_prova)

sizeFactors(spe_prova) <- spe_prova$Area_um / median(spe_prova$Area_um)
spe_prova <- logNormCounts(spe_prova, name="normalized_by_area")
assayNames(spe_prova)

df <- data.frame(colData(spe_prova), spatialCoords(spe_prova))
p1 <- ggplot(df, aes(x = labels_singler, y = library_size)) + 
  ggtitle("Library size factors") 

p2 <- ggplot(df, aes(x = labels_singler, y = sizeFactor)) + 
  ggtitle("Area-derived factors") 

(p1 + p2) &
  geom_boxplot(aes(fill = labels_singler)) & # No need to specify `x` again
  labs(x = "Cell type", y = "Size factor") &
  scale_fill_manual(values = unname(pals::tableau20())) &
  theme_bw() & 
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# zoom on y
(p1 + p2) &
  geom_boxplot(aes(fill = labels_singler)) &
  labs(x = "Cell type", y = "Size factor") &
  scale_fill_manual(values = unname(pals::tableau20())) &
  
  # ADD THIS LINE to set the y-axis limits without removing data
  # It sets the upper limit to 5 and leaves the lower limit to be determined automatically.
  coord_cartesian(ylim = c(NA, 5)) &
  
  theme_bw() & 
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# spatial variation ----

# most penalized genes by differences between the 2 normalization methods
mean_libs <- rowMeans(assays(spe_prova)$logcounts)
mean_area <- rowMeans(assays(spe_prova)$normalized_by_area)
diff <- abs(mean_libs - mean_area)
ord <- order(diff, decreasing=TRUE)
top_diff <- diff[ord[1:10]]
names(top_diff)

# vettori conteggio dei 2 geni con maggiori differenze tra le due normalizzazioni
mx_ls <- as.matrix(t(logcounts(spe_prova)[ord[1:5], , drop = FALSE]))
mx_area <- assay(spe_prova, "normalized_by_area")[ord[1:5], ] |> 
  t() |> as.matrix()

colnames(mx_ls) <- paste0(colnames(mx_ls), "_ls")
colnames(mx_area) <- paste0(colnames(mx_area), "_area")


colData(spe_prova) <- cbind(colData(spe_prova), mx_ls, mx_area)

p1 <- plotCentroids(spe_prova, colourBy = "Mt1_ls", size = 1) +
  theme(
    legend.key.width  = unit(0.5, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  #guides(color = guide_legend(override.aes = list(size = 5)))  +
  scale_y_reverse() +
  theme(
    axis.title = element_blank(), # Removes both x and y axis titles
    axis.text = element_blank(),  # Removes both x and y axis text labels
    axis.ticks = element_blank()  # Removes both x and y axis tick marks
  ) + labs(title="Library size normalization")
p2 <- plotCentroids(spe_prova, colourBy = "Mt1_area", size = 1) +
  theme(
    legend.key.width  = unit(0.5, "lines"),
    legend.key.height = unit(1, "lines")
  ) +
  #guides(color = guide_legend(override.aes = list(size = 5)))  +
  scale_y_reverse() +
  theme(
    axis.title = element_blank(), # Removes both x and y axis titles
    axis.text = element_blank(),  # Removes both x and y axis text labels
    axis.ticks = element_blank()  # Removes both x and y axis tick marks
  ) + labs(title="Area normalization") 

(p1 / p2) &
  coord_equal() & 
  theme_void()


# Dimensionality Reduction ----

library(Banksy)
set.seed(20000229)

spe <- spe_list$blocco2_c26
rm(spe_list)

dec <- modelGeneVar(spe)
hvg <- getTopHVGs(dec, n=3e3)

spe_hvg <- spe[hvg, ]
spe_hvg

# non spatial
spe_hvg <- logNormCounts(spe_hvg)
spe_hvg <- runPCA(spe_hvg, ncomponents=20)

# spatially aware
spe_hvg <- computeBanksy(spe_hvg, assay_name="logcounts", k_geom=30, coord_names = c("x_um, y_um"), )
spe_hvg <- runBanksyPCA(spe_hvg, npcs=20, lambda=1)

reducedDimNames(spe_hvg) <- c("PCA_tx", "PCA_sp")

# umap on both
spe_hvg <- runUMAP(spe_hvg, dimred="PCA_tx", name="UMAP_tx")
spe_hvg <- runUMAP(spe_hvg, dimred="PCA_sp", name="UMAP_sp")

# PCA on spatial coords
df <- data.frame(colData(spe_hvg), spatialCoords(spe_hvg), 
                 reducedDim(spe_hvg, "PCA_M0_lam0.8"))
p2 <- lapply(paste0("PC", seq(3)), \(.) {
  fd <- df[order(abs(df[[.]])), ]
  ggplot(fd, aes(x, y, col=.data[[.]]))
}) |>
  wrap_plots(nrow=1) &
  geom_point(shape=16, stroke=0, size=0.5) &
  scale_color_gradientn(colors=pals::jet()) &
  coord_equal() & theme_void() & theme(
    legend.position="bottom",
    legend.key.width=unit(0.8, "lines"),
    legend.key.height=unit(0.4, "lines"))

p2

# umap comparison
.plot_umap <- \(dr) plotReducedDim(
  object=spe_hvg, dimred=dr, colour_by="labels_singler", 
  point_shape=16, point_size=0) + ggtitle(dr)

.plot_umap("UMAP_tx") +
  .plot_umap("UMAP_sp") +
  plot_layout(guides="collect") &
  guides(col=guide_legend(override.aes=list(alpha=1, size=2))) &
  theme_void() & theme(aspect.ratio=1, 
                       legend.key.size=unit(0, "lines"),
                       plot.title=element_text(hjust=0.5))

# clustering
g <- buildSNNGraph(spe_hvg, use.dimred="PCA_tx", type="jaccard", k=15, num.threads = 10)
# cluster using Leiden community detection algorithm
k <- igraph::cluster_leiden(g, objective_function="modularity", resolution=0.5)
table(spe_hvg$cluster_leiden <- factor(k$membership))

plotCentroids(spe_hvg, colourBy = "cluster_leiden", size = 0.7)

g <- buildSNNGraph(spe_hvg, use.dimred="PCA_sp", type="jaccard", k=15, num.threads = 10)
k <- igraph::cluster_leiden(g, objective_function="modularity", resolution=0.5)
table(spe_hvg$cluster_banksy <- factor(k$membership))

plotCentroids(spe_hvg, colourBy = "cluster_banksy", size = 0.7)


# Joint counts

pcs <- reducedDim(spe_hvg, "PCA")[, 1:10]
params <- bluster::KNNGraphParam(
  k=20,
  cluster.fun="leiden",
  cluster.args=list(
    resolution=0.3,
    objective_function="modularity"))
colData(spe_hvg)$cluster <- bluster::clusterRows(pcs, BLUSPARAM=params)

plotCentroids(spe_hvg, colourBy = "cluster", size = 0.7)


resJc <- spdep::joincount.multi(as.factor(spe_hvg$cluster), colGraph(spe_hvg, "knn6"))
resJc <- resJc[order(resJc[, "z-value"], decreasing=TRUE), ]
head(resJc, 20)