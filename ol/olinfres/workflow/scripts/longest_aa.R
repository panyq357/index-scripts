library(Biostrings)

config <- list(
    raw_aa = "rawdata/OL_proteins.fa",
    longest_aa = "results/longest_aa.fa"
)

raw_aa <- readAAStringSet(config$raw_aa)

names(raw_aa) <- names(raw_aa) |> sub("OL([^T]+)T([^.]+)\\..*", "OL\\1G\\2", x=_, perl=T)

aa <- raw_aa[grep("^OL", x=names(raw_aa))]

longest_aa <- split(aa, names(aa)) |>
    lapply(function(x) {x[order(width(x), decreasing=T)][1]}) |>
    unname() |> 
    do.call(c, args=_)

writeXStringSet(longest_aa, config$longest_aa)

