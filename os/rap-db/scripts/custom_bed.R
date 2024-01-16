library(GenomicRanges)

config <- list(
    gtf = snakemake@input$gtf,
    tx = snakemake@output$transcript
)

main <- function() {

    gtf <- rtracklayer::import(config$gtf)
    tx <- gtf[gtf$type == "transcript",]

    write.table(gr_to_bed_df(tx), file=config$tx, quote=F, sep="\t", row.names=F, col.names=F)

}

gr_to_bed_df <- function(gr) {
    df <- data.frame(
        seqnames=seqnames(gr),
        starts=start(gr)-1,
        ends=end(gr),
        names=rep(".", length(gr)),
        scores=rep(".", length(gr)),
        strands=strand(gr)
    )
    return(df)
}

main()
