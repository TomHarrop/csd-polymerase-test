#!/usr/bin/env bash

set -eux

singularity exec \
    --nv \
    guppy_3.4.1.sif \
    guppy_basecaller \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --qscore_filtering \
    --input_path data/longamp \
    --recursive \
    --save_path 01_basecalled/longamp \
    --barcode_kits "EXP-PBC096" \
    --num_barcode_threads 4 \
    --detect_mid_strand_barcodes \
    --device "cuda:0"

singularity exec \
    --nv \
    guppy_3.4.1.sif \
    guppy_basecaller \
    --flowcell FLO-MIN106 \
    --kit SQK-LSK109 \
    --qscore_filtering \
    --input_path data/q5 \
    --recursive \
    --save_path 01_basecalled/q5 \
    --barcode_kits "EXP-PBC096" \
    --num_barcode_threads 4 \
    --detect_mid_strand_barcodes \
    --device "cuda:0"


singularity exec \
    porechop_0.2.4.sif \
    porechop \
    -i 01_basecalled/longamp/pass \
    -o 02_porechop/longamp.fastq \
    --verbosity 1 \
    --threads 8 \
    --check_reads 1000 \
    --discard_middle

singularity exec \
    porechop_0.2.4.sif \
    porechop \
    -i 01_basecalled/q5/pass \
    -o 02_porechop/q5.fastq \
    --verbosity 1 \
    --threads 8 \
    --check_reads 1000 \
    --discard_middle


singularity exec \
    ngmlr_8d76779 \
    ngmlr \
    -r data/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
    -q 02_porechop/longamp.fastq \
    -o 03_ngmlr/longamp.sam \
    --rg-sm "longamp" \
    --rg-id "longamp" \
    -t 8 \
    -x ont

singularity exec \
    ngmlr_8d76779 \
    ngmlr \
    -r data/GCF_003254395.2_Amel_HAv3.1_genomic.fna \
    -q 02_porechop/q5.fastq \
    -o 03_ngmlr/q5.sam \
    --rg-sm "q5" \
    --rg-id "q5" \
    -t 8 \
    -x ont


samtools sort \
    -@ 8 \
    03_ngmlr/longamp.sam \
    > 04_bam/longamp.bam

samtools index \
    04_bam/longamp.bam

samtools sort \
    -@ 8 \
    03_ngmlr/q5.sam \
    > 04_bam/q5.bam

samtools index \
    04_bam/q5.bam