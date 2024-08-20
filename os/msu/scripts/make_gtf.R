library(rtracklayer)
library(GenomicRanges)

config <- list(
    gff = snakemake@input$gff,
    out_gtf = snakemake@output$out_gtf,
    out_gtf_rap_chr = snakemake@output$out_gtf_rap_chr
)

if (!dir.exists(dirname(config$out_gtf))) dir.create(dirname(config$out_gtf))

gff <- rtracklayer::import(config$gff)

transcript <- subset(gff, type == "mRNA")
transcript$gene_id <- transcript$Parent
transcript$transcript_id <- transcript$ID
transcript$type <- "transcript"
transcript <- subset(transcript, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

CDS <- subset(gff, type == "CDS")
CDS$transcript_id <- CDS$Parent
CDS$gene_id <- sub("LOC_Os(\\d{2})g(\\d{5}).*", "LOC_Os\\1g\\2", CDS$transcript_id)
CDS <- subset(CDS, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

five_prime_utr <- subset(gff, type == "five_prime_UTR")
five_prime_utr$transcript_id <- five_prime_utr$Parent
five_prime_utr$gene_id <- sub("LOC_Os(\\d{2})g(\\d{5}).*", "LOC_Os\\1g\\2", five_prime_utr$transcript_id)
five_prime_utr$type <- "five_prime_utr"
five_prime_utr <- subset(five_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

three_prime_utr <- subset(gff, type == "three_prime_UTR")
three_prime_utr$transcript_id <- three_prime_utr$Parent
three_prime_utr$gene_id <- sub("LOC_Os(\\d{2})g(\\d{5}).*", "LOC_Os\\1g\\2", three_prime_utr$transcript_id)
three_prime_utr$type <- "three_prime_utr"
three_prime_utr <- subset(three_prime_utr, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

exon <- subset(gff, type == "exon")
exon$transcript_id <- exon$Parent
exon$gene_id <- sub("LOC_Os(\\d{2})g(\\d{5}).*", "LOC_Os\\1g\\2", exon$transcript_id)
exon <- subset(exon, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

gene <- subset(gff, type == "gene")
gene$transcript_id <- NA
gene$gene_id <- gene$ID
gene <- subset(gene, select = c("source", "type", "score", "phase", "gene_id", "transcript_id"))

new_gtf <- c(transcript, CDS, five_prime_utr, three_prime_utr, exon, gene)

new_gtf <- subset(new_gtf, ! seqnames(new_gtf) %in% c("ChrSy", "ChrUn"))

seqlevels(new_gtf) <- sprintf("Chr%d", 1:12)

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

new_gtf_rap_chr <- new_gtf
seqlevels(new_gtf_rap_chr) <- seqlevels(new_gtf_rap_chr) |> sub("Chr", "", x=_) |> as.integer() |> sprintf("chr%02d", ...=_)

export(new_gtf, config$out_gtf, format="gtf")
export(new_gtf_rap_chr, config$out_gtf_rap_chr, format="gtf")

