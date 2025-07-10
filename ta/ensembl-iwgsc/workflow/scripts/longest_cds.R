library(Biostrings)

cds_all <- readDNAStringSet(snakemake@params$url)

# Remove N containing CDS.
cds_all <- cds_all[vcountPattern("N", cds_all) == 0]

names(cds_all) <- stringr::str_extract(names(cds_all), "Traes[^.]+")

order_idx <- sapply(cds_all, length) |> split(names(cds_all)) |> lapply(order, decreasing=T) |> unlist()

longest_cds <- cds_all[order_idx == 1]

writeXStringSet(longest_cds, snakemake@output[[1]])
