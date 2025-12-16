library(arrow)

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
