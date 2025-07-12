library(GenomicRanges)

unzip(
  zipfile = snakemake@input[[1]],
  files = "iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3",
  exdir = snakemake@output$temp_dir
)

gr <- rtracklayer::import(file.path(snakemake@output$temp_dir, "iwgsc_refseqv2.1_gene_annotation_200916", "iwgsc_refseqv2.1_annotation_200916_HC.gff3"))

gr$Parent <- as.character(gr$Parent)  # For unknown reason, rtracklayer::import recognize Parent as CharacterList, need to convert to character to prevent error.

mcols(gr)$gene_id <- NA

parent_to_id <- match(as.character(mcols(gr)$Parent), mcols(gr)$ID)

mcols(gr)$gene_id <- mcols(gr)$Parent[parent_to_id]
mcols(gr)$gene_id[gr$type == "gene"] <- mcols(gr)$ID[gr$type == "gene"]
mcols(gr)$gene_id[gr$type == "mRNA"] <- mcols(gr)$Parent[gr$type == "mRNA"]

mcols(gr)$transcript_id <- mcols(gr)$Parent
mcols(gr)$transcript_id[gr$type == "mRNA"] <- as.character(mcols(gr)$ID)[gr$type == "mRNA"]

mcols(gr)$type <- as.character(mcols(gr)$type)
mcols(gr)$type[gr$type == "mRNA"] <- "transcript"
mcols(gr)$type[gr$type == "three_prime_UTR"] <- "three_prime_utr"
mcols(gr)$type[gr$type == "five_prime_UTR"] <- "five_prime_utr"

is_gene_coding <- gr$gene_id %in% gr$gene_id[gr$type == "CDS"]
gr$gene_biotype <- "ncRNA"
gr$gene_biotype[is_gene_coding] <- "protein_coding"
is_tx_coding <- (gr$transcript_id %in% gr$transcript_id[gr$type == "CDS"]) & (gr$type != "gene")
gr$transcript_biotype <- NA
gr$transcript_biotype[gr$type != "gene"] <- "ncRNA"
gr$transcript_biotype[is_tx_coding] <- "protein_coding"

mcols(gr) <- mcols(gr)[c("source", "type", "score", "phase", "gene_id", "transcript_id", "gene_biotype", "transcript_biotype")]

rtracklayer::export(gr, snakemake@output$gtf, format="gtf")
