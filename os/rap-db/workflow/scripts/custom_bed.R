library(rtracklayer)
library(GenomicRanges)

gtf <- rtracklayer::import(snakemake@input$gtf)
fa <- Biostrings::readDNAStringSet(snakemake@input$fa)
seqlengths(gtf) <- sapply(fa[names(seqlengths(gtf))], length)

transcript <- subset(gtf, type == "transcript")
mcols(transcript)$name <- transcript$transcript_id
mcols(transcript) <- mcols(transcript)["name"]
export(transcript, snakemake@output$transcript, format="bed")

gene <- subset(gtf, type == "gene")
mcols(gene)$name <- gene$gene_id
mcols(gene) <- mcols(gene)["name"]
export(gene, snakemake@output$gene, format="bed")

exon <- subset(gtf, type == "exon")
# Assign exon number
mcols(exon)$exon_number <- split(exon, exon$transcript_id) |>
  lapply(function(tx_exon) {
    cat(sprintf("\rAssign exon number for %s", tx_exon$transcript_id[1]))
    if (all(strand(tx_exon) == "+")) {
      return(order(start(tx_exon)))
    } else if (all(strand(tx_exon) == "-")) {
      return(order(end(tx_exon), decreasing = TRUE))
    } else {
      return(NA)
    }
  }) |>
  unlist()
cat("\n")
mcols(exon)$name <- with(mcols(exon), sprintf("%s.exon%d", transcript_id, exon_number))
mcols(exon) <- mcols(exon)["name"]
export(exon, snakemake@output$exon, format="bed")

upstream_3kb <- subset(gtf, type == "transcript") |>
  promoters(upstream=3000, downstream=0) |>
  trim()
mcols(upstream_3kb)$name <- upstream_3kb$transcript_id
mcols(upstream_3kb) <- mcols(upstream_3kb)["name"]
export(upstream_3kb, snakemake@output$upstream_3kb, format="bed")
