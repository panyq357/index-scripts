raw <- Biostrings::readDNAStringSet(snakemake@input[[1]])

genome <- raw[1:12]

gwh_ids <- sub("(\\S+)\\s.*", "\\1", names(genome))
ori_seq_ids <- sub(".*OriSeqID=(\\S+)\\s.*", "\\1", names(genome))

names(genome) <- ori_seq_ids

chr_order <- sub("Chr", "", ori_seq_ids) |> as.integer() |> order()
genome <- genome[chr_order]

Biostrings::writeXStringSet(genome, snakemake@output$genome)

chr_id_converter <- data.frame(gwh_ids, ori_seq_ids)

readr::write_tsv(chr_id_converter, snakemake@output$chr_id_converter)
