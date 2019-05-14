#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import gzip
import argparse
import gc

def handle_fastq(fq1_file, fq2_file, set1):
    # setU: 全集, set2: 储存重复的ID, seq_dict: 序列为key, ID列表为值
    with gzip.open(fq1_file) as f1, gzip.open(fq2_file) as f2:
        fq1 = f1.readlines()
        fq2 = f2.readlines()
    with open('stat.txt', 'w') as fout:
        fout.write('TPR\tFPR\tPrecision\tF1\tRange\tLength\n')
        for length in xrange(6, 151, 2):
            for pos in xrange(1, 146, 2):
                setU, set2, seq_dict = set(), set(), {}
                if pos+length-1 > 150:
                    break
                for n, line in enumerate(fq1):
                    if n % 4 == 0:
                        ID = line.split()[0][1:]
                        setU.add(ID)
                    if n % 4 == 1:
                        seq = line[pos-1:pos+length-1] + \
                                fq2[n][pos-1:pos+length-1]
                        if seq in seq_dict:
                            seq_dict[seq].append(ID)
                        else:
                            seq_dict[seq] = [ID]
                for seq, id_list in seq_dict.iteritems():
                    if len(id_list) > 1:
                        set2 |= set(id_list)
                # 计算各指标
                TP = len(set2 & set1)
                FN = len(set1 - set2)
                FP = len(set2 - set1)
                TN = len(setU - (set2 | set1))
                TPR = TP/float(TP+FN)
                FPR = FP/float(FP+TN)
                Precision = TP/float(TP+FP)
                F1 = 2*Precision*TPR/(Precision+TPR)
                rang = '{0}-{1}'.format(pos, pos+length-1)
                fout.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                    TPR, FPR, Precision, F1, rang, length*2))
                fout.flush()
                # 释放内存
                del set2
                del seq_dict
                gc.collect()

def handle_bam(bam_file):
    dup_set, set1 = set(), set()
    # 只考虑FLAG 1187, 1107, 1171, 1123为duplicates
    dupflags = {'1187','1107','1171','1123'}
    # 第一次读文件，记录被markduplicates的reads对的位点和方向
    with os.popen('samtools view {}'.format(bam_file)) as f:
        for line in f:
            ls = line.split()
            # ctg和pos确定一端read的位置，tlen确定另一端的位置和方向
            (seqid, flag, ctg, pos), tlen = ls[:4], ls[8]
            if flag in dupflags and tlen > 0:
                dup_set.add((ctg, pos, tlen))
    # 第二次读取文件，只要reads对位点和方向与markduplicates
    # 的reads对相同的都统计到set1
    with os.popen('samtools view {}'.format(bam_file)) as f:
        for line in f:
            ls = line.split()
            (seqid, flag, ctg, pos), tlen = ls[:4], ls[8]
            if (ctg, pos, tlen) in dup_set:
                set1.add(seqid)
    return set1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bam_file',
            help='picard MarkDuplicates BAM file')
    parser.add_argument('read1',
            help='read1 FASTQ file')
    parser.add_argument('read2',
            help='read2 FASTQ file')
    args = parser.parse_args()

    set1 = handle_bam(args.bam_file)
    handle_fastq(args.read1, args.read2, set1)

if __name__ == '__main__':
    main()

