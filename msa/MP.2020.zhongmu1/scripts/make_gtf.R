library(rtracklayer)
library(GenomicRanges)

config <- list(
    gff = snakemake@input[[1]],
    out_gtf = snakemake@output[[1]]
)

if (!dir.exists(dirname(config$out_gtf))) dir.create(dirname(config$out_gtf))

raw_gff <- rtracklayer::import(config$gff)

raw_gff <- subset(raw_gff, seqnames(raw_gff) %in% sprintf("Chr%s", 1:8))

transcript <- subset(raw_gff, type == "mRNA")
transcript$gene_id <- transcript$Parent
transcript$transcript_id <- transcript$ID
transcript$type <- "transcript"
transcript <- subset(transcript, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

CDS <- subset(raw_gff, type == "CDS")
CDS$transcript_id <- CDS$Parent
CDS$gene_id <- sub("(MsG.*)\\.T.*", "\\1", CDS$transcript_id)
CDS <- subset(CDS, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

five_prime_utr <- subset(raw_gff, type == "five_prime_UTR")
five_prime_utr$transcript_id <- five_prime_utr$Parent
five_prime_utr$gene_id <- sub("(MsG.*)\\.T.*", "\\1", five_prime_utr$transcript_id)
five_prime_utr$type <- "five_prime_utr"
five_prime_utr <- subset(five_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

three_prime_utr <- subset(raw_gff, type == "three_prime_UTR")
three_prime_utr$transcript_id <- three_prime_utr$Parent
three_prime_utr$gene_id <- sub("(MsG.*)\\.T.*", "\\1", three_prime_utr$transcript_id)
three_prime_utr$type <- "three_prime_utr"
three_prime_utr <- subset(three_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

exon <- subset(raw_gff, type == "exon")
exon$transcript_id <- exon$Parent
exon$gene_id <- sub("(MsG.*)\\.T.*", "\\1", exon$transcript_id)
exon$type <- "exon"
exon <- subset(exon, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

tx_by_gene <- split(transcript, as.character(transcript$gene_id))
gene <- unlist(reduce(tx_by_gene))
gene$gene_id <- names(gene)
names(gene) <- NULL
gene$transcript_id <- NA
gene$source <- "calculated"
gene$type <- "gene"
gene$phase <- NA
gene$score <- NA

new_gtf <- c(transcript, CDS, five_prime_utr, three_prime_utr, exon, gene)

new_gtf <- sort(new_gtf, by = ~ width, decreasing=T)
new_gtf <- sort(new_gtf, by = ~ start, decreasing=F)
new_gtf <- sort(new_gtf, by = ~ seqnames, decreasing=F)

is_gene_coding <- !is.na(match(as.character(new_gtf$gene_id), as.character(CDS$gene_id)))
new_gtf$gene_biotype <- "ncRNA"
new_gtf$gene_biotype[is_gene_coding] <- "protein_coding"
is_tx_coding <- !is.na(match(as.character(new_gtf$transcript_id), as.character(CDS$transcript_id)))
new_gtf$transcript_biotype <- "ncRNA"
new_gtf$transcript_biotype[is_tx_coding] <- "protein_coding"
new_gtf$transcript_biotype[new_gtf$type == "gene"] <- NA

export(new_gtf, config$out_gtf, format="gtf")

