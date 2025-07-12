
rule download_genome:
    output:
        "rawdata/iwgsc_refseqv2.1_assembly.fa.zip"
    params:
        url = "https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/iwgsc_refseqv2.1_assembly.fa.zip"
    shell:
        "mkdir -p rawdata; cd rawdata; aria2c -c {params.url}"


rule download_gene_annotation:
    output:
        "rawdata/iwgsc_refseqv2.1_gene_annotation_200916.zip"
    params:
        url = "https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseqv2.1_gene_annotation_200916.zip"
    shell:
        "mkdir -p rawdata; cd rawdata; aria2c -c {params.url}"



rule download_functional_annotation:
    output:
        "rawdata/iwgsc_refseqv2.1_functional_annotation.zip"
    params:
        url = "https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseqv2.1_functional_annotation.zip"
    shell:
        "mkdir -p rawdata; cd rawdata; aria2c -c {params.url}"
