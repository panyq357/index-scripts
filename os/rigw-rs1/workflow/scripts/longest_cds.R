library(Biostrings)

config <- list(
    raw_cds = snakemake@input$raw_cds,
    longest_cds = snakemake@output$longest_cds
)

raw_cds <- readDNAStringSet(config$raw_cds)

names(raw_cds) <- names(raw_cds) |> sub("([^t]+)t([^\\-]+)-.*", "\\1g\\2", x=_, perl=T)

longest_cds <- split(raw_cds, names(raw_cds)) |>
    lapply(function(x) {
        x[which(order(x) == 1)]
        return(x[1])
    }) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_cds, config$longest_cds)

