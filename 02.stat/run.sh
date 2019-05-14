#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2019-05-13 16:06

ln -s ../01.markdup/fq1.filter.fastq.gz
ln -s ../01.markdup/fq2.filter.fastq.gz
ln -s ../01.markdup/out.markdup.bam

python result_stat.py out.markdup.bam fq1.filter.fastq.gz fq2.filter.fastq.gz
