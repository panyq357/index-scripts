from pathlib import Path


rule download_genome:
    output:
        "rawdata/IRGSP-1.0_genome.fasta.gz"
    params:
        url = config["download_links"]["genome"]
    shell:
        "cd rawdata ; wget {params.url}"


gff_url = config["download_links"]["gff"]

rule download_gff:
    output:
        locus = "rawdata/IRGSP-1.0_representative/locus.gff",
        transcripts = "rawdata/IRGSP-1.0_representative/transcripts.gff",
        transcripts_exon = "rawdata/IRGSP-1.0_representative/transcripts_exon.gff"
    params:
        url = gff_url,
        file = Path(gff_url).name
    shell:
        "cd rawdata ; wget {params.url} ; tar -xzf {params.file}"


rule download_kasalath:
    output:
        "rawdata/kasalath_genome/all.fasta"
    params:
        url = config["download_links"]["kasalath"]
    shell:
        "cd rawdata ; wget {params.url} ; tar -xzf kasalath_genome.tar.gz"


rule download_annotation:
    output:
        Path("rawdata") / Path(config["download_links"]["annotation"]).name
    params:
        url = config["download_links"]["annotation"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_id_converter:
    output:
        Path("rawdata") / Path(config["download_links"]["id_converter"]).name
    params:
        url = config["download_links"]["id_converter"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_protein:
    output:
        Path("rawdata") / Path(config["download_links"]["protein"]).name
    params:
        url = config["download_links"]["protein"]
    shell:
        "cd rawdata ; wget {params.url}"


rule download_transcript:
    output:
        Path("rawdata") / Path(config["download_links"]["transcript"]).name
    params:
        url = config["download_links"]["transcript"]
    shell:
        "cd rawdata ; wget {params.url}"

