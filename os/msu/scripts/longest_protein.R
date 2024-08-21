library(Biostrings)

config <- list(
    raw_protein = snakemake@input$raw_protein,
    longest_protein = snakemake@output$longest_protein
)

raw_protein <- readAAStringSet(config$raw_protein)

names(raw_protein) <- names(raw_protein) |> sub("([^.]+)\\..*", "\\1", x=_, perl=T)

protein <- raw_protein[grep("^LOC", names(raw_protein)),]

longest_protein <- split(protein, names(protein)) |>
    lapply(function(x) {x[order(width(x), decreasing=T)][1]}) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_protein, config$longest_protein)

