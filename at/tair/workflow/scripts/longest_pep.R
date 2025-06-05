fa <- Biostrings::readAAStringSet(snakemake@input[[1]])
names(fa) <- sub("(AT\\dG\\d{5})\\..*", "\\1", names(fa))
Biostrings::writeXStringSet(fa, snakemake@output[[1]])
