library(purrr)
library(ggplot2)
library(SpatialExperiment)
#library(scuttle)
#library(SpaceTrooper)

# negate a statement ----
`%nin%` <- negate(`%in%`)

# dot plot ----
dot_plot_paper <- function(df) {
  df <- df[df$pct_expr_combat > 0, ]  # Filter out zero values
  p <- ggplot(df, aes(x = gene, y = celltype)) +
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
  return(p)
}

# marker dot plot

dot_plot_marker <- function(df, gene_to_plot) {
  
  # 1. Filter for the specific gene passed to the function
  df_filtered <- df[df$gene == gene_to_plot, ]
  
  # 2. Filter out rows where the expression percentage is zero (as preferred)
  df_filtered <- df_filtered[df_filtered$pct_expr > 0, ]
  
  # 3. Build the plot with the new aesthetics and labels
  p <- ggplot(df_filtered, aes(x = celltype, y = condition)) +
    geom_point(aes(size = pct_expr, color = avg_expr)) +
    scale_color_gradient(low = "grey", high = "red") +
    labs(
      title = gene_to_plot, 
      y = "Condition", 
      x = "Cell type",
      color = "Avg expression", 
      size = "% expr >0"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# mapping labels from sublist to complete list ----
mappingLabels <- function(spe_list, ann_list) {
  # Your existing loop is perfect
  for (name in names(spe_list)) {
    if (name %in% names(ann_list)) {
      
      spe <- spe_list[[name]]
      ann_df <- ann_list[[name]]
      
      match_indices <- match(colnames(spe), ann_df$cell_id)
      spe$labels_rctd <- ann_df$first_type[match_indices]
      spe_list[[name]] <- spe
      
      cat("Mapped 'labels_rctd' for object:", name, "\n")
      
    } else {
      cat("Skipping object:", name, "(no corresponding annotation found)\n")
    }
  }

  return(spe_list)
}

# reading polygons from parquet file to sf objects
read_shapes <- function(poly_path, samples, shape_name = "nuclei_boundaries"){
  poly_list <- map(samples, function(sample){
    poly <- st_read_parquet(paste0(poly_path, sample, ".zarr/shapes/", shape_name,"/shapes.parquet"))
    names(poly)[names(poly) == "X__index_level_0__"] <- "cell_id"
    rownames(poly) <- poly$cell_id
    poly
  })
}

# adding polygons to spe_obj
add2sfe <- function(spe_obj, polygons, MARGIN = NULL, name = "filtered_nuclei"){
  # filtering nuclei in the sample and reordering
  target_ids <- colData(spe_obj)[["cell_id"]]
  
  # Subset poly_list$blocco1_c26STAT3 to only rows with cell_id in target_ids
  keep <- polygons$cell_id %in% target_ids
  poly_filtered <- polygons[keep, ]
  
  # Reorder so cell_id matches order in target_ids
  ord <- match(target_ids, poly_filtered$cell_id)
  poly_filtered <- poly_filtered[ord, ]
  
  # (Optional) assign rownames to match colnames of spe_obj
  rownames(poly_filtered) <- colnames(spe_obj)
  
  # (Optional) check alignment
  stopifnot(all(poly_filtered$cell_id == target_ids))
  
  # add to the sfe
  # MARGIN: 1 = righe (geni)
  #         2 = colonne (nuclei)
  dimGeometry(spe_obj, name,  MARGIN = MARGIN) <- poly_filtered
  return(spe_obj)
}


# Function to clean the repeated coldata 
clean_coldata <- function(spe) {
  # Ensure the object has colData
  stopifnot(is(object = spe, "SpatialExperiment"))
  
  # Extract and process colData
  col_data <- colData(spe)
  col_names <- colnames(col_data)
  
  # Find duplicates (columns with the same name)
  unique_names <- unique(col_names)
  to_keep <- rep(TRUE, ncol(col_data)) # Keeps track of which columns to retain
  
  # Iterate over unique column names
  for (name in unique_names) {
    indices <- which(col_names == name)
    
    if (length(indices) > 1) {
      # Case 1: Columns are identical
      identical_cols <- sapply(indices, function(i) {
        all(col_data[[i]] == col_data[[indices[1]]])
      })
      
      if (all(identical_cols)) {
        # Keep only the first instance if all are identical
        to_keep[indices[-1]] <- FALSE
      } else {
        # Case 2: Columns with the same name but different content
        # Keep the one with the highest index (drop lower indices)
        to_keep[indices[-length(indices)]] <- FALSE
      }
    }
  }
  
  # Subset colData to only keep the valid columns
  clean_col_data <- col_data[, to_keep, drop = FALSE]
  
  # Update the original SpatialExperiment object with cleaned colData
  colData(spe) <- clean_col_data
  return(spe)
}