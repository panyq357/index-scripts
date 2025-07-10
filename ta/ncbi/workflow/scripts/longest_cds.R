cds_all <- Biostrings::readDNAStringSet(snakemake@params$url)

gene_names <- stringr::str_extract(names(cds_all), "(?<=\\[gene=)[^\\]]+")

order_idx <- unname(cds_all) |> sapply(length) |> split(gene_names) |> lapply(order, decreasing=TRUE) |> unlist()

out <- cds_all[order_idx == 1] |> setNames(gene_names[order_idx == 1])

Biostrings::writeXStringSet(out, snakemake@output[[1]])
