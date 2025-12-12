df <- readr::read_tsv(snakemake@input[[1]])[c("Locus_ID", "Oryzabase Gene Symbol Synonym(s)", "Description")]
df <- dplyr::distinct(df)

df <- df[match(unique(df[[1]]), df[[1]]), ]

readr::write_tsv(df, snakemake@output[[1]], col_names=FALSE)
