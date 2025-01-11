
rule download_build4_genome:
    output:
        expand("rawdata/OsGenome_RAP2/chr{num}.seq.masked", num=["01","02","03","04","06","07","08","09","10","11","12"])
    params:
        url = "https://rapdb.dna.affrc.go.jp/download/archive/build4/OsGenome_RAP2.tar.gz"
    shell:
        '''
        mkdir -p rawdata ; cd rawdata
        wget {params.url}
        tar -xf OsGenome_RAP2.tar.gz
        '''


rule download_build4_gff:
    output:
        "rawdata/gff_RAP2/rep.gff"
    params:
        url = "https://rapdb.dna.affrc.go.jp/download/archive/build4/gff_RAP2.tar.gz"
    shell:
        '''
        mkdir -p rawdata ; cd rawdata
        wget {params.url}
        tar -xf gff_RAP2
        '''


rule make_build4_genome:
    input:
        expand("rawdata/OsGenome_RAP2/chr{num}.seq.masked", num=["01","02","03","04","06","07","08","09","10","11","12"])
    output:
        fa = "results/genome/os.build4.genome.fa"
    shell:
        '''
        cat {input} > {output.fa}
        '''

rule make_build4_gene_bed:
    input:
        "rawdata/gff_RAP2/rep.gff"
    output:
        "results/custom_bed/build4_gene.bed"
    script:
        "../scripts/make_build4_gene_bed.R"
