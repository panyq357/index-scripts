library(GenomicRanges)


gff <- rtracklayer::import(snakemake@input$gff)
chr_id_converter <- readr::read_tsv(snakemake@input$chr_id_converter) |>
  with(setNames(ori_seq_ids, gwh_ids))


# Filter out features not in Chr1-12.
gff <- gff[seqnames(gff) %in% names(chr_id_converter), ]


# Rename chromosome names.
seqlevels(gff) <- seqlevels(gff)[seqlevels(gff) %in% names(chr_id_converter)]
seqlevels(gff) <- chr_id_converter[seqlevels(gff)]
seqlevels(gff) <- sprintf("Chr%d", 1:12)  # Reorder levels for sortting.


gene <- gff[gff$type == "gene", ]
gene$gene_id <- as.character(gene$ID)
mcols(gene) <- mcols(gene)[c("source", "type", "score", "phase", "gene_id")]

transcript <- subset(gff, type == "mRNA")
transcript$transcript_id <- as.character(transcript$ID)
transcript$gene_id <- as.character(transcript$Parent)
mcols(transcript) <- mcols(transcript)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]

sub_tx_feature_names <- setdiff(levels(gff$type), c("gene", "mRNA"))
sub_tx_feature_list <- sub_tx_feature_names |>
  lapply(function(feature_name) {
    feature <- subset(gff, type == feature_name)
    feature$transcript_id <- as.character(feature$Parent)
    feature$gene_id <- transcript$gene_id[match(feature$transcript_id, transcript$transcript_id)]
    mcols(feature) <- mcols(feature)[c("source", "type", "score", "phase", "transcript_id", "gene_id")]
    return(feature)
  }) |>
  setNames(sub_tx_feature_names) |>
  GRangesList()


new_gtf <- unlist(c(GRangesList(list(gene, transcript)), sub_tx_feature_list))
names(new_gtf) <- NULL

# Sort new_gtf
new_gtf <- new_gtf[order(new_gtf$type, decreasing=TRUE), ]
new_gtf <- new_gtf[order(width(new_gtf), decreasing=TRUE), ]
new_gtf <- new_gtf[order(start(new_gtf), decreasing=FALSE), ]
new_gtf <- new_gtf[order(seqnames(new_gtf), decreasing=FALSE), ]


levels(new_gtf$type) <- c("gene", "transcript", "five_prime_utr", "exon", "CDS", "three_prime_utr")


# Add biotype based on have CDS or not.
new_gtf$gene_biotype <- "ncRNA"
new_gtf$gene_biotype[new_gtf$gene_id %in% sub_tx_feature_list$CDS$gene_id] <- "protein_coding"
new_gtf$transcript_biotype <- NA
new_gtf$transcript_biotype[new_gtf$type != "gene"] <- "ncRNA"
new_gtf$transcript_biotype[new_gtf$type != "gene" & new_gtf$transcript_id %in% sub_tx_feature_list$CDS$transcript_id] <- "protein_coding"


rtracklayer::export(new_gtf, snakemake@output$gtf, format="gtf")
