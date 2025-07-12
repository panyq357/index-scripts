readr::read_csv(snakemake@input[[1]])[c("ZS11 Gene ID", "AtGI/Name", "Description")] |>
  readr::write_tsv(snakemake@output[[1]], col_names = FALSE)


