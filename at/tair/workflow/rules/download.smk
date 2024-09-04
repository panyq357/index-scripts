
rule download_genome:
    output:
        "rawdata/TAIR10_chr_all.fas.gz"
    params:
        url = "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz"
    shell:
        "cd rawdata ; aria2c {params.url}"


rule download_gff:
    output:
        "rawdata/TAIR10_GFF3_genes.gff"
    params:
        url = "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"
    shell:
        "cd rawdata ; aria2c {params.url}"


rule download_description:
    output:
        "rawdata/TAIR10_functional_descriptions"
    params:
        url = "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/TAIR10_genome_release/TAIR10_functional_descriptions"
    shell:
        "cd rawdata ; aria2c {params.url}"


rule download_pep:
    output:
        "rawdata/TAIR10_pep_20110103_representative_gene_model"
    params:
        url = "https://www.arabidopsis.org/api/download-files/download?filePath=Proteins/TAIR10_protein_lists/TAIR10_pep_20110103_representative_gene_model"
    shell:
        "cd rawdata ; aria2c {params.url}"

