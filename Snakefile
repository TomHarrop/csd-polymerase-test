#!/usr/bin/env python3

import multiprocessing

# GLOBALS

raw_reads = 'data/csd-polymerase-test.tar.gz'

threads = multiprocessing.cpu_count()

guppy = 'shub://TomHarrop/ont-containers:guppy_3.4.1'
ngmlr = 'shub://TomHarrop/align-utils:ngmlr_8d76779'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'

rule target:
    input:
        expand('output/040_ngmlr/{pol}.bam.bai',
               pol=['longamp', 'q5'])

rule index:
    input:
        'output/040_ngmlr/{pol}.bam'
    output:
        'output/040_ngmlr/{pol}.bam.bai'
    singularity:
        samtools
    shell:
        'samtools index {input}'

rule sort:
    input:
        'output/040_ngmlr/{pol}.sam'
    output:
        'output/040_ngmlr/{pol}.bam'
    log:
        'output/logs/sort_{pol}.log'
    singularity:
        samtools
    shell:
        'samtools sort '
        '{input} '
        '> {output} '
        '2> {log}'


rule ngmlr:
    input:
        reads = 'output/030_porechop/{pol}.fastq',
        ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna',
        fai = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna.fai'
    output:
        pipe('output/040_ngmlr/{pol}.sam')
    log:
        'output/logs/ngmlr_{pol}.log'
    threads:
        threads - 1
    singularity:
        ngmlr
    shell:
        'ngmlr '
        '-r {input.ref} '
        '-q {input.reads} '
        '--rg-sm {wildcards.pol} '
        '--rg-id {wildcards.pol} '
        '-t {threads} '
        '-x ont '
        '>> {output} '
        '2> {log}'

rule porechop:
    input:
        'output/020_fastq/{pol}/pass'
    output:
        'output/030_porechop/{pol}.fastq'
    log:
        'logs/porechop_{pol}.log'
    threads:
        threads // 2
    singularity:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--check_reads 1000 '
        '--discard_middle '
        '&> {log}'

rule basecall:
    input:
        'output/010_reads/{pol}'
    output:
        directory('output/020_fastq/{pol}/pass')
    log:
        'output/logs/guppy_{pol}.log'
    params:
        wd = 'output/020_fastq/{pol}'
    singularity:
        guppy
    threads:
        threads // 2
    shell:
        'guppy_basecaller '
        '--flowcell FLO-MIN106 '
        '--kit SQK-LSK109 '
        '--qscore_filtering '
        '--input_path {input} '
        '--recursive '
        '--save_path {params.wd} '
        '--barcode_kits "EXP-PBC096" '
        '--num_barcode_threads {threads} '
        '--detect_mid_strand_barcodes '
        '--device auto '
        '&> {log}'


rule unzip:
    input:
        raw_reads
    output:
        temp(directory('output/010_reads/longamp')),
        temp(directory('output/010_reads/q5'))
    params:
        wd = 'output/010_reads'
    singularity:
        samtools
    shell:
        'tar -zxf {input} '
        '-C {params.wd} '
        '--strip-components 1'