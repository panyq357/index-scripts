library(Biostrings)

config <- list(
    raw_cds = snakemake@input$raw_cds,
    longest_cds = snakemake@output$longest_cds
)

raw_cds <- readDNAStringSet(config$raw_cds)

names(raw_cds) <- names(raw_cds) |> sub(".*\\[gene=([^\\]]+)\\].*", "\\1", x=_, perl=T)

cds <- raw_cds[grep("^LOC", x=names(raw_cds))]

longest_cds <- split(cds, names(cds)) |>
    lapply(function(x) {x[order(width(x), decreasing=T)][1]}) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_cds, config$longest_cds)

