raw_cds <- Biostrings::readDNAStringSet(snakemake@input$cds)


# Remove CDS not in Chr1-12.
ori_seq_ids <- sub(".*OriSeqID=(\\S+)", "\\1", names(raw_cds))
cds <- raw_cds[grep("Chr\\d+", ori_seq_ids), ]


names(cds) <- sub(".*OriGeneID=(\\S+)\\s.*", "\\1", names(cds))


cds <- cds[order(Biostrings::width(cds), decreasing=TRUE)]
cds <- cds[!duplicated(names(cds))]


Biostrings::writeXStringSet(cds, snakemake@output$longest_cds)
