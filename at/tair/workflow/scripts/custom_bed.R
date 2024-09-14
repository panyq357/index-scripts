library(rtracklayer)
library(GenomicRanges)

config <- list(
    gtf = snakemake@input$gtf,

    transcript = snakemake@output$transcript,
    gene = snakemake@output$gene
)

main <- function() {

    gtf <- import(config$gtf)

    transcript <- gtf[gtf$type == "transcript",]
    mcols(transcript)$name <- transcript$transcript_id
    mcols(transcript) <- mcols(transcript)["name"]  # Keep only name in mcols, since rtracklayer can't export other columns.

    gene <- gtf[gtf$type == "gene",]
    mcols(gene)$name <- gene$gene_id
    mcols(gene) <- mcols(gene)["name"]

    export(transcript, config$transcript, format="bed")
    export(gene, config$gene, format="bed")

}

main()
