df <- readr::read_csv(snakemake@input[[1]])
df$symbol <- NA
df[df$f.type == "Human readable description", c("g2.identifier", "symbol", "f.name")] |> readr::write_tsv(snakemake@output[[1]], col_names=FALSE)

