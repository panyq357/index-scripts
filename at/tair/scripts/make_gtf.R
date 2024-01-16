library(rtracklayer)
library(GenomicRanges)

config <- list(
    raw_gtf = snakemake@input$raw_gtf,
    new_gtf = snakemake@output$new_gtf
)

gtf <- rtracklayer::import(config$raw_gtf)

# Convert factor to character.
mcols(gtf)$type <- as.character(mcols(gtf)$type)

# Rename type.
mcols(gtf)$type[mcols(gtf)$type == "five_prime_UTR"] <- "five_prime_utr"
mcols(gtf)$type[mcols(gtf)$type == "three_prime_UTR"] <- "three_prime_utr"
mcols(gtf)$type[mcols(gtf)$type == "mRNA"] <- "transcript"

# Add gene_biotype, "protein_coding" for genes have CDS, "ncRNA" for others.
gene_have_cds <- unique(gtf$gene_id[gtf$type == "CDS"])
is_gene_coding <- gtf$gene_id %in% gene_have_cds
gtf$gene_biotype <- "ncRNA"
gtf$gene_biotype[is_gene_coding] <- "protein_coding"

# Add transcript_biotype, "protein_coding" for transcripts have CDS, "ncRNA" for others.
tx_have_cds <- unique(gtf$transcript_id[gtf$type == "CDS"])
is_tx_coding <- gtf$transcript_id %in% tx_have_cds
gtf$transcript_biotype <- "ncRNA"
gtf$transcript_biotype[is_tx_coding] <- "protein_coding"

# Remove transcript_id and transcript_biotype of genes.
gtf$transcript_id[gtf$type == "gene"] <- NA
gtf$transcript_biotype[gtf$type == "gene"] <- NA

# Sort
gtf <- gtf[order(seqnames(gtf), start(gtf), end(gtf)),]

export(gtf, config$new_gtf, format="gtf")

