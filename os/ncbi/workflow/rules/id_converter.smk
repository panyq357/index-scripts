
rule longest_cds:
    input:
        raw_cds = "rawdata/GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna.gz"
    output:
        longest_cds = "results/longest_cds.fa"
    script:
        "../scripts/longest_cds.R"

rule diamond_blastx:
    input:
        longest_cds = "results/longest_cds.fa",
        d = "../rap-db/results/os.rap-db.dmnd"
    output:
        o = "results/diamond/ncbi.to.rap-db.tsv",
    shell:
        '''
        diamond blastx \
            -d {input.d} \
            -q {input.longest_cds} \
            -o {output.o} \
            --very-sensitive
        '''

rule id_converter:
    input:
        diamond_output = "results/diamond/ncbi.to.rap-db.tsv",
    output:
        best_hit = "results/id_converter/ncbi.to.rap-db.best-hit.tsv",
    script:
        "../scripts/id_converter.R"

