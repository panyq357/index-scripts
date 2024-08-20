library(Biostrings)

config <- list(
    raw_protein = snakemake@input$raw_protein,
    longest_protein = snakemake@output$longest_protein
)

raw_protein <- readAAStringSet(config$raw_protein)

names(raw_protein) <- names(raw_protein) |> sub("Os([^t]+)t([^-]+)-.*", "Os\\1g\\2", x=_, perl=T)

longest_protein <- split(raw_protein, names(raw_protein)) |>
    lapply(function(x) {
        x[which(order(x) == 1)]
        return(x[1])
    }) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_protein, config$longest_protein)

