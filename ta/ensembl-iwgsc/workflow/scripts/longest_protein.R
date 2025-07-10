library(Biostrings)

readDNAStringSet(snakemake@input[[1]]) |>
  translate() |>
  writeXStringSet(snakemake@output[[1]])
