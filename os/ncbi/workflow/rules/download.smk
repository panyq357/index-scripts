
rule download_cds:
    output:
        "rawdata/GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna.gz"
    params:
        outdir = "rawdata",
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/034/140/825/GCF_034140825.1_ASM3414082v1/GCF_034140825.1_ASM3414082v1_cds_from_genomic.fna.gz"
    shell:
        '''
        mkdir -p {params.outdir}
        cd {params.outdir}
        wget {params.url}
        '''
