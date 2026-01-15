
gff <- rtracklayer::import("rawdata/Brassica_napus.ZS11.v0.gene.gff3.gz")

# Sort GFF
gff <- sort(gff, by = ~ end, decreasing = FALSE)
gff <- sort(gff, by = ~ seqnames + start + type)

rtracklayer::export(gff, "results/make_gff/bn.bnir.zs11.v0.gene.gff", format="gff3")
