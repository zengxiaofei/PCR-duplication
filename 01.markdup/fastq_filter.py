#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import gzip

def handle_sam(sam_file):
    idset = set()
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            idset.add(line.split()[0])
    return idset

def handle_fastq(fastq_file, idset, out_file):
    with gzip.open(fastq_file) as fin, gzip.open(out_file, 'w') as fout:
        line = fin.readline()
        while line:
            if line.split()[0][1:] in idset:
                fout.write(line)
                fout.write(fin.readline())
                fout.write(fin.readline())
                fout.write(fin.readline())
            else:
                fin.readline()
                fin.readline()
                fin.readline()
            line = fin.readline()

def main():
    idset = handle_sam(sys.argv[1])
    handle_fastq(sys.argv[2], idset, sys.argv[4])
    handle_fastq(sys.argv[3], idset, sys.argv[5])

if __name__ == '__main__':
    main()
