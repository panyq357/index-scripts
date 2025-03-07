from pathlib import Path

# For each gene, keep its longest protein.
rule longest_protein:
    input:
        raw_protein = Path("rawdata") / Path(config["download_links"]["protein"]).name
    output:
        longest_protein = "results/longest_protein.fa"
    script:
        "../scripts/longest_protein.R"


# Make a diamond db, so that other assembly can use it to make id conveter.
rule diamond_makedb:
    input:
        "results/longest_protein.fa"
    output:
        "results/os.rap-db.dmnd"
    params:
        d = "results/os.rap-db"
    shell:
        "diamond makedb --in {input} -d {params.d}"


# This id converter is just a clean up of id converter from RAP-DB.
rule id_converter:
    input:
        Path("rawdata") / Path(config["download_links"]["id_converter"]).name
    output:
        "results/id_converter.txt"
    script:
        "../scripts/id_converter.R"
