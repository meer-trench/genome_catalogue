#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 11:02:18 2022

run_fastp.py -i sample.path.tsv -t {threads} -o sample.merged.fa.gz -f sample.unmerged.1.fa.gz -r sample.unmegred.2.fa.gz

@author: songz
"""

import argparse
import os
import sys
import hashlib
import datetime
# from zetaSeq import io as seqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Path to the sample path tav file.')
parser.add_argument('-w', '--raw', default='/dssg/home/acct-trench/trench/DATA/MEER_RAW/', help='Path to the raw data')
parser.add_argument('-t', '--threads', type=int, help='Number of threads.')
parser.add_argument('-o', '--output', help='Path to merged file.')
parser.add_argument('-f', '--forward', help='Path to unmerged r1 file.')
parser.add_argument('-r', '--reverse', help='Path to unmerged r2 file.')
parser.add_argument('-s', '--sequencer', choices=['mgi', 'other'], help='Sequencers')
parser.add_argument('-not_gz_out', action='store_true', help='Indicator for not compress output files.')
parser.add_argument('-fastq_out', action='store_true', help='Indicator for output fastq files instead of fasta.') 
parser.add_argument('-log', default='log.tsv', help='Path to log file.')
parser.add_argument('-d', '--dry_run', action='store_true', help='Indicator of dry-run.')
args=parser.parse_args()
input_file = args.input
raw_path = args.raw
if not raw_path.endswith('/'): raw_path += '/'
threads = args.threads
merged_file = args.output
unmerged_file = (args.forward, args.reverse)
mgi = False
if args.sequencer == 'mgi':
    mgi = True
    ad1 = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'
    ad2 = 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'
not_gz_out = args.not_gz_out
fastq_out = args.fastq_out
log_file = args.log
dry_run = args.dry_run

merge = False
if merged_file: merge = True


# Check if input file ends with .gz
if not not_gz_out:
    if merged_file:
        if merged_file.endswith('.gz'):
            pass
        else:
            print('Input file has to be .gz')
            sys.exit()
    if unmerged_file[0].endswith('.gz'):
        pass
    else:
        print('Input file has to be .gz')
        sys.exit()
    if unmerged_file[1].endswith('.gz'):
        pass
    else:
        print('Input file has to be .gz')
        sys.exit()


# Get FASTQ path
path= {}
with open(input_file, 'r') as f: # Each line is a path to a sequencing lane
    for line in f:
        line = line.strip('\n')
        content = []
        for file in os.listdir(raw_path + line):
            if file.endswith('.fq.gz') or file.endswith('.fastq'):
                content.append(file)
        content.sort()
        if len(content) != 2:
            print('Wrong file number.')
            sys.exit()
        path[line] = tuple(content)


# Run fastp
tmp = []
for key, value in path.items():
    print('Run fastp on path {0}'.format(key))
    hash_object = hashlib.md5((input_file+key+unmerged_file[0]).encode())
    tmp_seq = 'TMP/tmp_' + hash_object.hexdigest()
    tmp.append(tmp_seq)
    cmd = ['fastp']
    cmd += ['-i ' + raw_path + key + '/' + value[0]]
    cmd += ['-I ' + raw_path + key + '/' + value[1]]
    if merge: cmd += ['-m --merged_out ' + tmp_seq + '.merged.fq.gz']
    cmd += ['--out1 ' + tmp_seq + '.unmerged.1.fq.gz']
    cmd += ['--out2 ' + tmp_seq + '.unmerged.2.fq.gz']
    if mgi: cmd += ['--adapter_sequence ' + ad1]
    if mgi: cmd += ['--adapter_sequence_r2 ' + ad2]
    cmd += ['--length_required 100']
    cmd += ['--cut_front --cut_right -W 4 -M 20']
    cmd += ['-w ' + str(threads)]
    cmd += ['-j ' + log_file]
    print(' '.join(cmd))
    if not dry_run:
        os.system(' '.join(cmd))
    print('Finished fastp on path {0}'.format(key))


# Concat all finished file, and remove the temp files
if dry_run:
    if merge: print('Write merged to {0}'.format(merged_file))
    print('Write unmerged to {0} {1}'.format(unmerged_file[0], unmerged_file[1]))
else:
    if fastq_out:
        seqtk_cmd = 'seqtk seq'
    else:
        seqtk_cmd = 'seqtk seq -A'
    if merge: cmd1 = [seqtk_cmd, tmp[0] + '.merged.fq.gz', '>', merged_file.strip('.gz')]
    cmd2 = [seqtk_cmd, tmp[0] + '.unmerged.1.fq.gz', '>', unmerged_file[0].strip('.gz')]
    cmd3 = [seqtk_cmd, tmp[0] + '.unmerged.2.fq.gz', '>', unmerged_file[1].strip('.gz')]
    if merge: 
        print(' '.join(cmd1))
        if not dry_run: os.system(' '.join(cmd1))
    print(' '.join(cmd2))
    if not dry_run: os.system(' '.join(cmd2))
    print(' '.join(cmd3))
    if not dry_run: os.system(' '.join(cmd3))
    if len(tmp) > 1: # If more than one path presents
        for item in tmp[1:]:
            if merge: cmd1 = [seqtk_cmd, item + '.merged.fq.gz', '>>', merged_file.strip('.gz')]
            cmd2 = [seqtk_cmd, item + '.unmerged.1.fq.gz', '>>', unmerged_file[0].strip('.gz')]
            cmd3 = [seqtk_cmd, item + '.unmerged.2.fq.gz', '>>', unmerged_file[1].strip('.gz')]
            if merge:
                print(' '.join(cmd1))
                if not dry_run: os.system(' '.join(cmd1))
            print(' '.join(cmd2))
            if not dry_run: os.system(' '.join(cmd2))
            print(' '.join(cmd3))
            if not dry_run: os.system(' '.join(cmd3))


# Pigz all FASTA files
if not_gz_out:
    pass
else:
    if merge:
        cmd = ['pigz']
        cmd += ['-p ' + str(threads)]
        cmd += [merged_file.strip('.gz')]
        print(' ' .join(cmd))
        if not dry_run: os.system(' '.join(cmd))

    cmd = ['pigz']
    cmd += ['-p ' + str(threads)]
    cmd += [unmerged_file[0].strip('.gz')]
    print(' '.join(cmd))
    if not dry_run: os.system(' '.join(cmd))

    cmd = ['pigz']
    cmd += ['-p ' + str(threads)]
    cmd += [unmerged_file[1].strip('.gz')]
    print(' '.join(cmd))
    if not dry_run: os.system(' '.join(cmd))


# Remove TMP files
if dry_run:
    pass
else:
    for item in tmp:
        if merge: os.remove(item + '.merged.fq.gz')
        os.remove(item + '.unmerged.1.fq.gz')
        os.remove(item + '.unmerged.2.fq.gz')
print('TMP files removed')
if merge:
    print('QCed files in {0} {1} {2}'.format(merged_file, unmerged_file[0], unmerged_file[1]))
else:
    print('QCed files in {0} {1}'.format(unmerged_file[0], unmerged_file[1]))
