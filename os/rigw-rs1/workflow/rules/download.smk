
rule download_genome:
    params:
        url = "http://rice.hzau.edu.cn/rice_rs1/download_ext/MH63RS1.LNNK00000000.fsa.tar.gz"
    output:
        "rawdata/MH63RS1.LNNK00000000.fsa"
    shell:
        '''
        mkdir -p rawdata
        cd rawdata
        wget {params.url}
        tar -xzf MH63RS1.LNNK00000000.fsa.tar.gz
        '''


# After downloaded gff, a manual correction in line 1064446 needs to be done.
# A "\t" before "LOC_Os05g34730.1" needs to be changed to " ".
rule download_gff:
    params:
        url = "http://rice.hzau.edu.cn/rice_rs1/download_ext/MH63_chr.gff.tar.gz"
    output:
        "rawdata/MH63_chr_v2.gff"
    shell:
        '''
        mkdir -p rawdata
        cd rawdata
        wget {params.url}
        tar -xzf MH63_chr.gff.tar.gz
        '''


rule download_cds:
    params:
        url = "http://rice.hzau.edu.cn/rice_rs1/download_ext/MH63.RS1.CDS.fa.tar.gz"
    output:
        "rawdata/MH63.RS1.CDS.fa"
    shell:
        '''
        mkdir -p rawdata
        cd rawdata
        wget {params.url}
        tar -xzf MH63.RS1.CDS.fa.tar.gz
        '''

