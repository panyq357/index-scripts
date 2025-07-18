library(GenomicRanges)

gr <- rtracklayer::import(snakemake@input[[1]])

## Remove duplicated mRNA rows.
## Can't use GenomicRanges::duplicated, it will ignore type.
is_duplicated <- duplicated(sprintf("%s:%s-%s(%s) %s %s %s", seqnames(gr), start(gr), end(gr), strand(gr), gr$type, gr$ID, gr$Parent))
gr <- gr[!is_duplicated]

## Sort all features.
seqlevels(gr) <- sprintf("chr%02d", 1:12)
gr <- sort(gr, by = ~ end, decreasing = FALSE)
gr <- sort(gr, by = ~ seqnames + start + type)

## Match gene_id.
mcols(gr)$gene_id <- NA
parent_to_id <- match(as.character(mcols(gr)$Parent), mcols(gr)$ID)  # Match exon, CDS, five/three_prime_utr's Parent to their transcript's ID.
mcols(gr)$gene_id <- as.character(mcols(gr)$Locus_id)[parent_to_id]  # Use their transcript's Locus_id as their gene_id.
mcols(gr)$gene_id[gr$type == "gene"] <- mcols(gr)$ID[gr$type == "gene"]
mcols(gr)$gene_id[gr$type == "mRNA"] <- mcols(gr)$Locus_id[gr$type == "mRNA"]

## Get transcript_id.
mcols(gr)$transcript_id <- as.character(mcols(gr)$Parent)
mcols(gr)$transcript_id[gr$type == "mRNA"] <- as.character(mcols(gr)$ID)[gr$type == "mRNA"]

## Rename type to match ensembl GTF format.
mcols(gr)$type <- as.character(mcols(gr)$type)
mcols(gr)$type[gr$type == "mRNA"] <- "transcript"
mcols(gr)$type[gr$type == "three_prime_UTR"] <- "three_prime_utr"
mcols(gr)$type[gr$type == "five_prime_UTR"] <- "five_prime_utr"

## Get gene_biotype and transcript_biotype for VEP annotation.
## Default is ncRNA.
## If a gene's gene_id matches any CDS's gene_id, it is a protein_coding gene.
## Same as transcript.
is_gene_coding <- gr$gene_id %in% gr$gene_id[gr$type == "CDS"]
gr$gene_biotype <- "ncRNA"
gr$gene_biotype[is_gene_coding] <- "protein_coding"
is_tx_coding <- (gr$transcript_id %in% gr$transcript_id[gr$type == "CDS"]) & (gr$type != "gene")
gr$transcript_biotype <- NA
gr$transcript_biotype[gr$type != "gene"] <- "ncRNA"
gr$transcript_biotype[is_tx_coding] <- "protein_coding"

## Remove annotation columns.
mcols(gr) <- mcols(gr)[c("source", "type", "score", "phase", "gene_id", "transcript_id", "gene_biotype", "transcript_biotype")]

rtracklayer::export(gr, snakemake@output[[1]], format="gtf")
