#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:57:34 2024

meer_snk_scripts/scripts/bwa_mem_100_samples.py {input.path} {output.mk} {threads}

bwa mem -t {threads} {params.index} {input.a1} {input.a2} | samtools sort -@ {threads} -o {output.bam} - 

He who loves to comment his code is unlikely to have bad luck.
@author: Song, Zewei
@contact: songzewei@genomics.cn
"""

import sys
import os

input_path = sys.argv[1]
output_mk = sys.argv[2]
threads = sys.argv[3]

sn = os.path.basename(input_path).split('.')[0]
index = 'data/contigs_for_binning/' + sn
sns = []
with open(input_path, 'r') as f:
    for line in f:
        line = line.strip('\n')
        sns.append(line)

for item in line:
    a1 = 'data/fastp_not_merged/' + item + '.unmerged_1.fastq'
    a2 = 'data/fastp_not_merged/' + item + '.unmerged_2.fastq'
    bam = 'data/alignments/' + item + '/' + item + '.sorted.bam'
    cmd = f'bwa mem -t {threads} {index} {a1} {a2} | samtools sort -@ {threads} -o {bam} -'
    print(cmd)
    os.system(cmd)