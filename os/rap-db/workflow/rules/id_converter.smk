
# For each gene, keep its longest protein.
rule longest_protein:
    input:
        raw_protein = "rawdata/IRGSP-1.0_protein_2024-07-12.fasta.gz"
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
        "rawdata/RAP-MSU_2023-09-07.txt.gz"
    output:
        "results/id_converter.txt"
    script:
        "scripts/id_converter.py"
