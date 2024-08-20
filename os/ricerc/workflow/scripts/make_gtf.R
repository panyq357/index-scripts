library(GenomicRanges)

config <- list(
    variety = snakemake@wildcards[["variety"]]
)

path_list <- list(
    tar = sprintf("rawdata/%s.IGDBv1.Allset.gff.tar.gz", config$variety),
    gff = sprintf("%s/%s.IGDBv1.Allset.gff", config$variety, config$variety),
    gtf = sprintf("results/gtf/%s.make.gtf", config$variety)
)

main <- function() {
    untar(path_list$tar)

    gff <- rtracklayer::import(path_list$gff)
    gff <- subset(gff, seqnames(gff) %in% sprintf("Chr%d", 1:12))

    gene <- subset(gff, type == "gene")
    gene$gene_id <- as.character(gene$ID)
    mcols(gene) <- mcols(gene)[c("source", "type", "score", "phase", "gene_id")]

    transcript <- subset(gff, type == "mRNA")
    transcript$transcript_id <- as.character(transcript$ID)
    transcript$gene_id <- as.character(transcript$Parent)
    transcript$type <- "transcript"
    mcols(transcript) <- mcols(transcript)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]

    sub_tx_feature_names <- c("exon", "CDS", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon")
    sub_tx_feature_list <- sub_tx_feature_names |>
        lapply(sub_tx_feature_pipeline, gff=gff, transcript=transcript) |>
        setNames(sub_tx_feature_names) |>
        GRangesList()

    # Sort new_gtf
    new_gtf <- unlist(c(GRangesList(list(gene, transcript)), sub_tx_feature_list))
    names(new_gtf) <- NULL
    new_gtf <- new_gtf[order(width(new_gtf), decreasing=T),]
    new_gtf <- new_gtf[order(start(new_gtf), decreasing=F),]
    new_gtf <- new_gtf[order(seqnames(new_gtf), decreasing=F),]

    # Add biotype based on have CDS or not.
    new_gtf$gene_biotype <- "ncRNA"
    new_gtf$gene_biotype[new_gtf$gene_id %in% sub_tx_feature_list$CDS$gene_id] <- "protein_coding"
    new_gtf$transcript_biotype <- NA
    new_gtf$transcript_biotype[new_gtf$type != "gene"] <- "ncRNA"
    new_gtf$transcript_biotype[new_gtf$type != "gene" & new_gtf$transcript_id %in% sub_tx_feature_list$CDS$transcript_id] <- "protein_coding"

    if (!dir.exists(dirname(path_list$gtf))) dir.create(dirname(path_list$gtf), recursive=T)
    rtracklayer::export(new_gtf, path_list$gtf, format="gtf")
}

sub_tx_feature_pipeline <- function(feature_name, gff, transcript) {
    feature <- subset(gff, type == feature_name)
    feature$transcript_id <- as.character(feature$Parent)
    feature$gene_id <- get_gene_id_by_transcript_id(feature$transcript_id, transcript)
    mcols(feature) <- mcols(feature)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]
    return(feature)
}

get_gene_id_by_transcript_id <- function(transcript_id, ref_gr) {
    gene_id <- ref_gr$gene_id[match(transcript_id, ref_gr$transcript_id)]
    return(gene_id)
}

main()
