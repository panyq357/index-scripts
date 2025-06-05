library(GenomicRanges)

# Rules for finding best hit
# - Lowest EValue
# - Identity > 80

df <- read.table(snakemake@input$diamond_output, col.names=c(
  "Query", "Target", "Identity", "Length",
  "Mismatch", "Gap", "QueryStart", "QueryEnd",
  "TargetStart", "TargetEnd", "EValue", "BitScore"
))

best_hit_pairs <- split(df, df$Query) |>
  sapply(function(df) {
#    df <- subset(df, Identity > 80)
    if (nrow(df) < 1) {
      return(NA)
    } else {
      return(df$Target[order(df$EValue, decreasing=F)][1])
    }
  })

write.table(data.frame(best_hit_pairs), snakemake@output$best_hit, row.names=T, col.names=F, sep="\t", quote=F)
