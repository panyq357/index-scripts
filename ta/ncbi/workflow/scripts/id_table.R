gtf <- rtracklayer::import(snakemake@params$url, format="gtf")

id_table <- GenomicRanges::mcols(gtf)[c("gene_id", "transcript_id", "protein_id")] |> as.data.frame() |> dplyr::distinct()

id_table[id_table == ""] <- NA

readr::write_tsv(id_table, snakemake@output[[1]])
