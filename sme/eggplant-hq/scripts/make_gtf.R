library(rtracklayer)
library(GenomicRanges)

config <- list(
    gff = "rawdata/02.SME-HQ.gff.gz",
    chromosomes = sprintf("E%02d", 1:12),
    out_gtf = "results/gtf/sm.eggplant-hq.make.gtf"
)


main <- function() {

    if (!dir.exists(dirname(config$out_gtf))) dir.create(dirname(config$out_gtf))

    gff <- rtracklayer::import(config$gff)

    transcript <- subset(gff, type == "mRNA")
    transcript$gene_id <- transcript$Parent
    transcript$transcript_id <- transcript$ID
    transcript$type <- "transcript"
    transcript <- subset(transcript, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

    CDS <- subset(gff, type == "CDS")
    CDS$transcript_id <- CDS$Parent
    CDS$gene_id <- get_gene_id_by_transcript_id(as.character(CDS$transcript_id), transcript)
    CDS <- subset(CDS, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

    five_prime_utr <- subset(gff, type == "five_prime_UTR")
    five_prime_utr$transcript_id <- five_prime_utr$Parent
    five_prime_utr$gene_id <- get_gene_id_by_transcript_id(as.character(five_prime_utr$transcript_id), transcript)
    five_prime_utr$type <- "five_prime_utr"
    five_prime_utr <- subset(five_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

    three_prime_utr <- subset(gff, type == "three_prime_UTR")
    three_prime_utr$transcript_id <- three_prime_utr$Parent
    three_prime_utr$gene_id <- get_gene_id_by_transcript_id(as.character(three_prime_utr$transcript_id), transcript)
    three_prime_utr$type <- "three_prime_utr"
    three_prime_utr <- subset(three_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

    exon <- subset(gff, type == "exon")
    exon$transcript_id <- exon$Parent
    exon$gene_id <- get_gene_id_by_transcript_id(as.character(exon$transcript_id), transcript)
    exon <- subset(exon, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

    new_gtf <- c(transcript, CDS, five_prime_utr, three_prime_utr, exon)
    new_gtf <- new_gtf[seqnames(new_gtf) %in% config$chromosomes,]

    new_gtf <- sort(new_gtf, by = ~ width, decreasing=T)
    new_gtf <- sort(new_gtf, by = ~ start, decreasing=F)
    new_gtf <- sort(new_gtf, by = ~ seqnames, decreasing=F)

    export(new_gtf, config$out_gtf, format="gtf")
}

#' Rationale:
#'   ref_gr contains ranges that have both gene_id and transcript_id,
#'   so use match and subset to get gene_id by transcript_id.
#' Note:
#'   This function is vectorized.
get_gene_id_by_transcript_id <- function(transcript_id, ref_gr) {
    gene_id <- ref_gr$gene_id[match(transcript_id, ref_gr$transcript_id)]
    return(gene_id)
}

main()
