readr::read_tsv(snakemake@input[[1]])[c("Locus_ID", "Oryzabase Gene Symbol Synonym(s)", "Description")] |>
  readr::write_tsv(snakemake@output[[1]], col_names=F)
