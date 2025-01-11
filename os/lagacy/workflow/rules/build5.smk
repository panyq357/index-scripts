rule download_build5_genome:
    output:
        "rawdata/IRGSPb5.fa.masked.gz"
    params:
        url = "https://rapdb.dna.affrc.go.jp/download/archive/build5/IRGSPb5.fa.masked.gz"
    shell:
        '''
        mkdir -p rawdata ; cd rawdata
        wget {params.url}
        '''


rule download_build5_gff:
    output:
        "rawdata/GFF3_representative/RAP_locus.gff3"
    params:
        url = "https://rapdb.dna.affrc.go.jp/download/archive/build5/GFF3_representative.tar.gz"
    shell:
        '''
        mkdir -p rawdata ; cd rawdata
        wget {params.url}
        tar -xf GFF3_representative.tar.gz
        '''


rule make_build5_genome:
    input:
        "rawdata/IRGSPb5.fa.masked.gz"
    output:
        fa = "results/genome/os.build5.genome.fa"
    shell:
        '''
        zcat {input} | sed -E 's/build05r1.fasta//g' > {output.fa}
        '''

rule make_build5_gene_bed:
    input:
        "rawdata/GFF3_representative/RAP_locus.gff3"
    output:
        "results/custom_bed/build5_gene.bed"
    script:
        "../scripts/make_build5_gene_bed.R"
