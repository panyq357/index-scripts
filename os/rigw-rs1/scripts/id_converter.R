library(GenomicRanges)

config <- list(
    ref_gtf = snakemake@input$ref_gtf,
    diamond_output = snakemake@input$diamond_output,
    variety_gtf = snakemake@input$variety_gtf,

    out = snakemake@output$best_hit,

    diamond_output_col_names = c(
        "Query", "Target", "Identity", "Length",
        "Mismatch", "Gap", "QueryStart", "QueryEnd",
        "TargetStart", "TargetEnd", "EValue", "BitScore"
    )
)

main <- function() {

    df <- read.table(config$diamond_output, col.names=config$diamond_output_col_names)
    ref_gtf <- rtracklayer::import(config$ref_gtf)
    variety_gtf <- rtracklayer::import(config$variety_gtf)

    best_hit_pairs <- split(df, df$Query) |>
        sapply(pick_1st_evalue)

    gene_to_transcript_order_by_width <- get_gene_to_transcript_order_by_width(variety_gtf)

    id_converter <- list()
    for (gene_id in names(gene_to_transcript_order_by_width)) {
        found <- F
        for (transcript_id in gene_to_transcript_order_by_width[[gene_id]]) {
            best_hit_transcript <- best_hit_pairs[transcript_id]
            if (length(best_hit_transcript) > 0) {
                id_converter[[gene_id]] <- best_hit_transcript
                found <- T
                break
            }
        }
        if (found) {
            next
        } else {
            id_converter[[gene_id]] <- NA
        }
    }
    id_converter <- data.frame(names(id_converter), get_gene_id_by_transcript_id(id_converter, ref_gtf))

    write.table(id_converter, config$out, row.names=F, col.names=F, sep="\t", quote=F)
}

pick_1st_evalue <- function(df) {
    return(df$Target[order(df$EValue, decreasing=F)[1]])
}

get_gene_id_by_transcript_id <- function(transcript_id, ref_gr) {
    gene_id <- ref_gr$gene_id[match(transcript_id, ref_gr$transcript_id)]
    return(gene_id)
}

get_gene_to_transcript_order_by_width <- function(gr) {
    transcript <- subset(gr, type == "transcript")
    by_gene_tx_ids <- split(transcript$transcript_id, transcript$gene_id)
    by_gene_tx_width_order <- split(width(transcript), transcript$gene_id) |>
        lapply(function(x) order(x, decreasing=T))
    transcript_id_order_by_width <- mapply(function(x, y) x[y], x=by_gene_tx_ids, y=by_gene_tx_width_order)
    return(transcript_id_order_by_width)
}

main()
