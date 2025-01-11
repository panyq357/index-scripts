library(GenomicRanges)

gff <- rtracklayer::import(snakemake@input[[1]])

gene <- subset(gff, type=="gene")

mcols(gene) <- data.frame(name=mcols(gene)[["ID"]])

rtracklayer::export(gene, snakemake@output[[1]])
