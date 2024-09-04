library(GenomicRanges)

config <- list(
    gff_dir = "rawdata/gff",
    out_gtf = "results/gtf/ol.make.gtf"
)

if (!dir.exists(dirname(config$out_gtf))) dir.create(dirname(config$out_gtf), recursive=T)

gff <- list.files(config$gff_dir, full.name=T) |>
    lapply(rtracklayer::import) |>
    GRangesList() |>
    unlist()

transcript <- subset(gff, type == "mRNA")
transcript$gene_id <- transcript$Parent
transcript$transcript_id <- transcript$ID
transcript$type <- "transcript"
transcript <- subset(transcript, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

tx_to_gene <- function(tx) sub("OL(\\d{2})T(\\d{6}).*", "OL\\1G\\2", tx)

CDS <- subset(gff, type == "CDS")
CDS$transcript_id <- CDS$Parent
CDS$gene_id <- tx_to_gene(CDS$transcript_id)
CDS <- subset(CDS, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

five_prime_utr <- subset(gff, type == "five_prime_UTR")
five_prime_utr$transcript_id <- five_prime_utr$Parent
five_prime_utr$gene_id <- tx_to_gene(five_prime_utr$transcript_id)
five_prime_utr$type <- "five_prime_utr"
five_prime_utr <- subset(five_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

three_prime_utr <- subset(gff, type == "three_prime_UTR")
three_prime_utr$transcript_id <- three_prime_utr$Parent
three_prime_utr$gene_id <- tx_to_gene(three_prime_utr$transcript_id)
three_prime_utr$type <- "three_prime_utr"
three_prime_utr <- subset(three_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

exon <- c(CDS, five_prime_utr, three_prime_utr)
exon$type <- "exon"
exon <- subset(exon, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

gene <- subset(gff, type == "gene")
gene$transcript_id <- NA
gene$gene_id <- gene$ID
gene <- subset(gene, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

new_gtf <- c(transcript, CDS, five_prime_utr, three_prime_utr, exon, gene)

seqlevels(new_gtf) <- sub("chromosome([^_]+)_.*", "\\1", seqlevels(new_gtf)) |> as.integer() |> as.character()

new_gtf$type <- factor(new_gtf$type, levels=c("gene", "transcript", "exon", "five_prime_utr", "CDS", "three_prime_utr"))

new_gtf <- sort(new_gtf, by = ~ type, decreasing=F)
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

rtracklayer::export(new_gtf, config$out_gtf, format="gtf")

