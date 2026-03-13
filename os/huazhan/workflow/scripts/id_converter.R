diamond_output_col_names <- c(
  "Query", "Target", "Identity", "Length",
  "Mismatch", "Gap", "QueryStart", "QueryEnd",
  "TargetStart", "TargetEnd", "EValue", "BitScore"
)

df <- read.table(snakemake@input$diamond_output, col.names=diamond_output_col_names)

df <- df[order(df$Identity, df$Length, decreasing=TRUE), ]
df <- df[!duplicated(df$Query), ]

df[c("Query", "Target", "Identity")] |>
  write.table(snakemake@output$best_hit, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
