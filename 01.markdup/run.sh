#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2019-05-09 23:25

ln -s ../00.data/GCA_002507455.1_ASM250745v1_genomic.fna
ln -s ../00.data/SRR8203406.sra
pigz -p 8 SRR8203406_1.fastq
pigz -p 8 SRR8203406_2.fastq

# mapping
fq1=SRR8203406_1.fastq.gz
fq2=SRR8203406_2.fastq.gz
ref=GCA_002507455.1_ASM250745v1_genomic.fna

bwa index ${ref}
bwa mem -t 8 ${ref} ${fq1} ${fq2} > out.sam
samtools view -@ 8 -bS out.sam -o out.bam

# selection of properly paired reads
samtools flagstat -@ 8 out.bam # properly paired: 98.96%
samtools view -@ 8 out.bam -h -f 83 > out.filter.sam
samtools view -@ 8 out.bam -f 163 >> out.filter.sam
samtools view -@ 8 out.bam -f 99 >> out.filter.sam
samtools view -@ 8 out.bam -f 147 >> out.filter.sam
samtools view -bS -@ 8 out.filter.sam -o out.filter.bam
samtools flagstat -@ 8 out.filter.bam # properly paired: 100%

# filtering FASTQ files
python fastq_filter.py \ 
    out.filter.sam \
    ${fq1} ${fq2} \
    fq1.filter.fastq.gz \
    fq2.filter.fastq.gz
rm *.sam

# mark duplicates
java -jar ~/software/picard/build/libs/picard.jar \
    I=out.filter.bam \
    O=out.markdup.bam \
    M=markdup.metrics

