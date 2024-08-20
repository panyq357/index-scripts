library(biomaRt)

ensembl_osativa <- useEnsemblGenomes(biomart = "plants_mart", dataset = "osativa_eg_gene")

write.csv(listFilters(ensembl_osativa), "temp.csv")

pathway_to_gene <- read.table("results/osa_pathway/pathway_to_gene.tsv", col.names=c("Pathway", "Gene"))

res <- getBM(
    attributes = c("entrezgene_accession", "ensembl_gene_id"),
    filters = "entrezgene_accession",
    values = unique(pathway_to_gene[["Gene"]]),
    mart = ensembl_osativa
)

res <- subset(res, duplicated(entrezgene_accession))
res[order(res[[1]]),]
