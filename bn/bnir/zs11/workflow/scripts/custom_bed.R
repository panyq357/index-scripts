library(rtracklayer)
library(GenomicRanges)

gtf <- import(snakemake@input$gtf)

transcript <- gtf[gtf$type == "transcript",]
mcols(transcript)$name <- transcript$transcript_id
mcols(transcript) <- mcols(transcript)["name"]  # Keep only name in mcols, since rtracklayer can't export other columns.

gene_range_list <- split(gtf, gtf$gene_id) |> range()
gene <- unlist(gene_range_list)
gene$name <- names(gene_range_list)

export(transcript, snakemake@output$transcript, format="bed")
export(gene, snakemake@output$gene, format="bed")
