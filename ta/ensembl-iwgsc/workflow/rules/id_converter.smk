rule longest_cds:
    output:
        "results/longest_cds.fasta"
    params:
        url = "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/triticum_aestivum/cds/Triticum_aestivum.IWGSC.cds.all.fa.gz"
    script:
        "scripts/longest_cds.R"


rule longest_protein:
    input:
        "results/longest_cds.fasta"
    output:
        "results/longest_protein.fasta"
    script:
        "scripts/longest_protein.R"


rule diamond_blastx_to_rapdb:
    input:
        longest_cds = "results/longest_cds.fasta",
        d = "../../os/rap-db/results/os.rap-db.dmnd"
    output:
        o = "results/diamond/ensembl-iwgsc.to.rap-db.tsv",
    shell:
        '''
        diamond blastx \
            -d {input.d} \
            -q {input.longest_cds} \
            -o {output.o} \
            --very-sensitive
        '''


rule id_converter_to_rapdb:
    input:
        diamond_output = "results/diamond/ensembl-iwgsc.to.rap-db.tsv",
    output:
        best_hit = "results/id_converter/ensembl-iwgsc.to.rap-db.best-hit.tsv",
    script:
        "scripts/id_converter.R"


rule diamond_makedb:
    input:
        "results/longest_protein.fasta"
    output:
        "results/ta.ensembl-iwgsc.dmnd"
    params:
        d = "results/ta.ensembl-iwgsc"
    shell:
        "diamond makedb --in {input} -d {params.d}"


rule diamond_blastx_to_msu:
    input:
        longest_cds = "results/longest_cds.fasta",
        d = "../../os/msu/results/os.msu.dmnd"
    output:
        o = "results/diamond/ensembl-iwgsc.to.msu.tsv",
    shell:
        '''
        diamond blastx \
            -d {input.d} \
            -q {input.longest_cds} \
            -o {output.o} \
            --very-sensitive
        '''


rule id_converter_to_msu:
    input:
        diamond_output = "results/diamond/ensembl-iwgsc.to.msu.tsv",
    output:
        best_hit = "results/id_converter/ensembl-iwgsc.to.msu.best-hit.tsv",
    script:
        "scripts/id_converter.R"

