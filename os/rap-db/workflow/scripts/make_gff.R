library(GenomicRanges)

gr <- rtracklayer::import(snakemake@input[[1]])

## Remove duplicated mRNA rows.
## Can't use GenomicRanges::duplicated, it will ignore type.
is_duplicated <- duplicated(sprintf("%s:%s-%s(%s) %s %s %s", seqnames(gr), start(gr), end(gr), strand(gr), gr$type, gr$ID, gr$Parent))
gr <- gr[!is_duplicated]

## Sort all features.
seqlevels(gr) <- sprintf("chr%02d", 1:12)
gr <- sort(gr, by = ~ end, decreasing = FALSE)
gr <- sort(gr, by = ~ seqnames + start + type)

gr$Parent[gr$type == "mRNA"] <- gr$Locus_id[gr$type == "mRNA"]

mcols(gr) <- mcols(gr)[c("source", "type", "score", "phase", "ID", "Parent", "Name")]


## Correct "mRNA" with no CDS to "ncRNA"
gr$type <- as.character(gr$type)
have_cds <- gr[gr$type == "mRNA"]$ID %in% as.character(gr[gr$type == "CDS"]$Parent)
gr[gr$type == "mRNA"][!have_cds]$type <- "ncRNA"

rtracklayer::export(gr, snakemake@output[[1]], format="gff3")
