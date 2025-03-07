df <- readr::read_tsv(snakemake@input[[1]], col_names=c("RAP", "MSU"))

df$MSU <- gsub("(LOC_Os\\d{2}g\\d{5})\\.\\d", "\\1", df$MSU)

df <- tidyr::separate_longer_delim(df, "MSU", ",") |>
  dplyr::distinct()

df$RAP[df$RAP == "None"] <- NA
df$MSU[df$MSU == "None"] <- NA

readr::write_tsv(df, snakemake@output[[1]], col_names=F)
