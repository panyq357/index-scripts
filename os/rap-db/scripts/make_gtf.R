library(rtracklayer)
library(GenomicRanges)

config <- list(
    gff = "rawdata/IRGSP-1.0_representative/transcripts.gff",
    out_gtf = "results/gtf/os.rap-db.make.gtf"
)

if (!dir.exists(dirname(config$out_gtf))) dir.create(dirname(config$out_gtf))

tx_gff <- rtracklayer::import(config$gff)

transcript <- subset(tx_gff, type == "mRNA")
transcript$gene_id <- transcript$Locus_id
transcript$transcript_id <- transcript$ID
transcript$type <- "transcript"
transcript <- subset(transcript, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

CDS <- subset(tx_gff, type == "CDS")
CDS$transcript_id <- CDS$Parent
CDS$gene_id <- sub("Os(\\d{2})t(\\d{7}).*", "Os\\1g\\2", CDS$transcript_id)
CDS <- subset(CDS, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

five_prime_utr <- subset(tx_gff, type == "five_prime_UTR")
five_prime_utr$transcript_id <- five_prime_utr$Parent
five_prime_utr$gene_id <- sub("Os(\\d{2})t(\\d{7}).*", "Os\\1g\\2", five_prime_utr$transcript_id)
five_prime_utr$type <- "five_prime_utr"
five_prime_utr <- subset(five_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

three_prime_utr <- subset(tx_gff, type == "three_prime_UTR")
three_prime_utr$transcript_id <- three_prime_utr$Parent
three_prime_utr$gene_id <- sub("Os(\\d{2})t(\\d{7}).*", "Os\\1g\\2", three_prime_utr$transcript_id)
three_prime_utr$type <- "three_prime_utr"
three_prime_utr <- subset(three_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

exon <- c(CDS, five_prime_utr, three_prime_utr)
exon_by_tx <- split(exon, as.character(exon$transcript_id))
exon_by_tx <- reduce(exon_by_tx)
exon <- unlist(exon_by_tx)
exon$transcript_id <- rep(names(exon_by_tx), elementNROWS(exon_by_tx))
exon$gene_id <- sub("Os(\\d{2})t(\\d{7}).*", "Os\\1g\\2", exon$transcript_id)
exon$source <- "calculated"
exon$type <- "exon"
exon$phase <- NA
exon$score <- NA
names(exon) <- NULL

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

seqlevels(new_gtf) <- sprintf("chr%02d", 1:12)
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

