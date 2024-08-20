
config <- list(
    pathway_to_gene = snakemake@input$pathway_to_gene,
    ncbi_to_rapdb = snakemake@input$ncbi_to_rapdb,
    pathway_to_rapdb = snakemake@output$pathway_to_rapdb
)

pathway_to_gene <- read.table(config$pathway_to_gene, col.names=c("Pathway", "Gene"))

ncbi_to_rapdb <- read.table(config$ncbi_to_rapdb, col.names=c("NCBI", "RAP"))
ncbi_to_rapdb <- with(ncbi_to_rapdb, setNames(RAP, NCBI))

pathway_to_gene$Gene <- sprintf("LOC%s", pathway_to_gene$Gene)
pathway_to_gene$RAP <- ncbi_to_rapdb[pathway_to_gene$Gene]

pathway_to_rapdb <- subset(pathway_to_gene, ! is.na(RAP), select=c("Pathway", "RAP"))

write.table(pathway_to_rapdb, config$pathway_to_rapdb, row.names=F, col.names=F, sep="\t", quote=F)

