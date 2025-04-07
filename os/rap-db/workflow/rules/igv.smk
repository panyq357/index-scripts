rule compress_genome:
    input:
        "results/genome/os.rap-db.genome.fa"
    output:
        gz = "results/igv/os.rap-db.genome.fa.gz",
        fai = "results/igv/os.rap-db.genome.fa.gz.fai"
    shell:
        '''
        bgzip < {input} > {output.gz}
        samtools faidx {output.gz}
        '''


rule sort_gff:
    input:
        "rawdata/IRGSP-1.0_representative/transcripts.gff"
    output:
        gff = "results/igv/os.rap-db.transcripts.sorted.gff.gz",
        tbi = "results/igv/os.rap-db.transcripts.sorted.gff.gz.tbi"
    shell:
        '''
        cat {input} \
        | grep -v '^#' \
        | sort -k1,1 -k4,4n -k5,5n -t$'\\t' \
        | bgzip > {output.gff}
        tabix {output.gff}
        '''


rule copy_json:
    input:
        fai = "results/igv/os.rap-db.genome.fa.gz.fai",
        tbi = "results/igv/os.rap-db.transcripts.sorted.gff.gz.tbi",
        json = "os.rap-db.genome.json"
    output:
        "results/igv/os.rap-db.genome.json"
    shell:
        "cp {input.json} {output}"
