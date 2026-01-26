library(SpatialExperiment)
library(arrow)
library(dplyr)
library(tibble)

# Estrazione annotazione tipi cellulari da R per python ----

# estraggo da filtered_sfe (senza zone infiammate quindi) l'indice e l'annotazione dei nuclei

colData(sfe_list$blocco2_c26)[c("cell_id", "labels_singler")]

df_list <- purrr::map(sfe_list, function(sfe) {
  filt <- sfe[, !sfe$to_discard]
  cd <- colData(filt)
  return(cd[c("cell_id", "labels_singler")])
})

head(df_list)

output_dir <- "/mnt/europa/valerio/data/data_tables/annotation_results_filtered/my_segmentation"

for (name in names(df_list)) {
  df <- as.data.frame(df_list[[name]])
  write_parquet(df, file.path(output_dir, paste0(name, ".parquet")))
}

# Spaceranger ----


head(c(colnames(sr_spe$blocco1_c26STAT3),colData(sr_spe$blocco1_c26STAT3)[["labels_singler"]]))
prova1 <- colnames(sr_spe$blocco1_c26STAT3)
prova2 <- colData(sr_spe$blocco1_c26STAT3)[["labels_singler"]]
prova3 <- cbind(prova1, prova2)
 
df_list <- purrr::map(sr_spe, function(spe) {
  cell_id <- colnames(spe)
  labels_singler <- colData(spe)[["labels_singler"]]
  cd <- as.data.frame(cbind(cell_id, labels_singler))
  return(cd)
})

output_dir <- "/mnt/europa/valerio/data/data_tables/annotation_results_filtered/spaceranger"

for (name in names(df_list)) {
  df <- as.data.frame(df_list[[name]])
  write_parquet(df, file.path(output_dir, paste0(name, ".parquet")))
}
# --------------------------------------------------------------------------

# Estrazione colonne informative per il filtraggio dei nuclei
# area_um e rilevamento outlier, log2countArea ecc

cols_extract_fun <- function(
    spe_big, 
    spe_small) {
  cols_to_save <- setdiff(names(as.data.frame(colData(spe_big))), names(as.data.frame(colData(spe_small))))
  if (length(cols_to_save) > 0) {
    
    # 3. Extract these columns
    # We usually want to keep row identifiers (like barcodes/sample_id) to map them back later
    # Assuming rownames are the barcodes
    df_to_save <- as.data.frame(colData(spe_big)) |>
      select(cell_id, all_of(cols_to_save))
    
    message("Dataframe with ", length(cols_to_save), " unique columns")
    return(df_to_save)
    
  } else {
    message("No unique columns found in spe_big relative to spe_small")
  }
}

output_dir <- "/mnt/europa/valerio/data/data_tables/spaceranger/nuclei_filtering_info"

walk2(spe_list_up, spe_list, function(spe_big, spe_small){
  df <- cols_extract_fun(spe_big, spe_small)
  name <- unique(spe_big$sample_id)
  write_parquet(df, file.path(output_dir, paste0(name, ".parquet")))
})

# my segmentation instead ----


output_dir <- "/mnt/europa/valerio/data/data_tables/my_segmentation/nuclei_filtering_info"

iwalk(spe_list, function(spe, name){
  df <- as.data.frame(colData(spe))
  write_parquet(df, file.path(output_dir, paste0(name, ".parquet")))
})
