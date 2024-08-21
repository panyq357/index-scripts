library(Biostrings)

config <- list(
    cds_tar = snakemake@input$cds_tar,
    longest_cds = snakemake@output$longest_cds,
    untar = snakemake@output$untar
)

untar(config$cds_tar)

raw_cds <- readDNAStringSet(config$untar)

names(raw_cds) <- names(raw_cds) |> sub("([^.]+\\.[^.]+)\\..*", "\\1", x=_, perl=T)

longest_cds <- split(raw_cds, names(raw_cds)) |>
    lapply(function(x) {x[order(width(x), decreasing=T)][1]}) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_cds, config$longest_cds)

