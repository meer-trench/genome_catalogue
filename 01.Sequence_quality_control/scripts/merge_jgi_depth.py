#!/bin/env python
import sys
import os

in_folder = sys.argv[1]
out_path = sys.argv[2]
if not in_folder.endswith('/'): in_folder += '/'

depth = {}
depth_file = [in_folder + f for f in os.listdir(in_folder)]
with open(depth_file[0], 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        if line[0] == 'contigName':
            depth[line[0]] = line[1:]
        else:
            depth[line[0]] = line[1:]
            depth[line[0]][0] = float(depth[line[0]][0])
            depth[line[0]][1] = float(depth[line[0]][1])
            depth[line[0]][2] = float(depth[line[0]][2])
print('Found {0} contigs'.format(len(depth)-1))

count = 2
for item in depth_file[1:]:
    print('Processing {0} ({1})'.format(item, count))
    with open(item, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if line[0] == 'contigName':
                depth[line[0]] += [line[3], line[4]]
            else:
                depth[line[0]][1] += float(line[2])
                depth[line[0]] += [float(line[3]), float(line[4])]
    count += 1

with open(out_path, 'w') as f:
    for key, value in depth.items():
        f.write('{0}\t{1}\n'.format(key, '\t'.join([str(i) for i in value])))
