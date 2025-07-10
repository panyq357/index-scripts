diamond_output_col_names = c(
  "Query", "Target", "Identity", "Length",
  "Mismatch", "Gap", "QueryStart", "QueryEnd",
  "TargetStart", "TargetEnd", "EValue", "BitScore"
)

df <- read.table(snakemake@input$diamond_output, col.names=diamond_output_col_names)

split_order <- with(df, split(Identity, Query)) |> lapply(order, decreasing=TRUE) |> unlist()

out <- df[split_order == 1, c("Query", "Target", "Identity")]

write.table(out, snakemake@output$best_hit, row.names=FALSE, col.names=F, sep="\t", quote=F)
