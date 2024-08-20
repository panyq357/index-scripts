library(GenomicRanges)

# Rules for finding best hit
# - Lowest EValue
# - Identity > 80

config <- list(
    diamond_output = snakemake@input$diamond_output,
    out = snakemake@output$best_hit,

    diamond_output_col_names = c(
        "Query", "Target", "Identity", "Length",
        "Mismatch", "Gap", "QueryStart", "QueryEnd",
        "TargetStart", "TargetEnd", "EValue", "BitScore"
    )
)


main <- function() {
    
    if (!dir.exists(dirname(config$out))) dir.create(dirname(config$out), recursive=T)

    df <- read.table(config$diamond_output, col.names=config$diamond_output_col_names)

    best_hit_pairs <- split(df, df$Query) |>
        sapply(function(df) {
            df <- subset(df, Identity > 80)
            if (nrow(df) < 1) {
                return(NA)
            } else {
                return(df$Target[order(df$EValue, decreasing=F)][1])
            }
        })

    write.table(data.frame(best_hit_pairs), config$out, row.names=T, col.names=F, sep="\t", quote=F)
}

main()
