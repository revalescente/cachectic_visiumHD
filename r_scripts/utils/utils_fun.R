library(purrr)
library(ggplot2)

# negate a statement ----
`%nin%` <- negate(`%in%`)

# dot plot ----
dot_plot_paper <- function(df)
{
  p1 <- ggplot(df, aes(x = gene, y = celltype)) +
    geom_point(aes(size = pct_expr_combat, color = avg_expr_combat)) +
    scale_color_gradient(low = "grey", high = "red") +
    theme_bw() +
    labs(
      x = "Gene",
      y = "Cell type",
      color = "Average expression",
      size = "% expressing"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid.major = element_line(color = "grey90")
    )
  return(p1)
}

# custom spatialQC metrics ----

vHD_spatialPerCellQC <- function(spe, micronConvFact=0.316, rmZeros=TRUE) {
  stopifnot(is(object=spe, "SpatialExperiment"))

  spe <- addPerCellQC(spe)

  if(!all(spatialCoordsNames(spe) %in% names(colData(spe)))) {
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
  }

  spnc <- spatialCoords(spe) * micronConvFact
  colnames(spnc) <-  paste0(spatialCoordsNames(spe), "_um")
  colData(spe) <- cbind.DataFrame(colData(spe), spnc)
  spe$Area_um <- spe$area * (micronConvFact^2)
  spe$AspectRatio <- spe$major_axis_length/spe$minor_axis_length
  if ("AspectRatio" %in% colnames(colData(spe))) {
    spe$log2AspectRatio <- log2(spe$AspectRatio) # not cosmx
  } else { warning("Missing aspect ratio in colData") }
  
  spe$CountArea <- spe$sum/spe$Area_um
  spe$log2CountArea <- log2(spe$CountArea)
  if (rmZeros) {
    if (sum(spe$sum==0) > 0) {
      message("Removing ", dim(spe[,spe$sum==0])[2],
              " cells with 0 counts!")
      spe <- spe[,!spe$sum==0]
    }
  }
  return(spe)
}



