rule faidx:
    input:
        "{prefix}.fa"
    output:
        "{prefix}.fa.fai"
    shell:
        "samtools faidx {input}"


rule CreateSequenceDictionary:
    input:
        "{prefix}.fa"
    output:
        "{prefix}.dict"
    shell:
        "gatk CreateSequenceDictionary -R {output.fa}"
