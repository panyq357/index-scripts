library(Biostrings)

raw_protein <- readAAStringSet(snakemake@input$raw_protein)

names(raw_protein) <- names(raw_protein) |> sub("Os([^t]+)t([^-]+)-.*", "Os\\1g\\2", x=_, perl=T)

longest_protein <- split(raw_protein, names(raw_protein)) |>
    lapply(function(x) {x[order(width(x), decreasing=T)][1]}) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_protein, snakemake@output$longest_protein)

