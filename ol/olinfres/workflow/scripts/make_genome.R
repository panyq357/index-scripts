library(Biostrings)

genome <- readDNAStringSet("rawdata/OL_genome.fa")

genome <- genome[grep("chromosome", names(genome))]

names(genome) <- sub("chromosome([^_]+)_.*", "\\1", names(genome)) |> as.integer()

if (!dir.exists("results/genome")) dir.create("results/genome", recursive=T)

writeXStringSet(genome, "results/genome/ol.make.genome.fa")
