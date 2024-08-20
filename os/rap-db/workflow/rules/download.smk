
rule download_genome:
    output:
        "rawdata/IRGSP-1.0_genome.fasta.gz"
    params:
        url = config["download_links"]["genome"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_gff:
    output:
        locus = "rawdata/IRGSP-1.0_representative/locus.gff",
        transcripts = "rawdata/IRGSP-1.0_representative/transcripts.gff",
        transcripts_exon = "rawdata/IRGSP-1.0_representative/transcripts_exon.gff"
    params:
        url = config["download_links"]["gff"]
    shell:
        "cd rawdata ; wget {params.url} ; tar -xzf IRGSP-1.0_representative_2023-09-07.tar.gz"


rule download_kasalath:
    output:
        "rawdata/kasalath_genome/all.fasta"
    params:
        url = config["download_links"]["kasalath"]
    shell:
        "cd rawdata ; wget {params.url} ; tar -xzf kasalath_genome.tar.gz"


rule download_annotation:
    output:
        "rawdata/IRGSP-1.0_representative_annotation_2023-09-07.tsv.gz"
    params:
        url = config["download_links"]["annotation"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_id_converter:
    output:
        "rawdata/RAP-MSU_2023-09-07.txt.gz"
    params:
        url = config["download_links"]["id_converter"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_protein:
    output:
        "rawdata/IRGSP-1.0_protein_2024-07-12.fasta.gz"
    params:
        outdir = "rawdata",
        url = "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_protein_2024-07-12.fasta.gz"
    shell:
        "mkdir -p {params.outdir} ; cd {params.outdir} ; wget {params.url}"

