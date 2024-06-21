#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 09:26:46 2024

meer_snk_scripts/scripts/random_draw_99_samples.py {input.path} {output.path}

He who loves to comment his code is unlikely to have bad luck.
@author: Song, Zewei
@contact: songzewei@genomics.cn
"""

import sys
import os
import pathlib
from random import sample

input_path = sys.argv[1]
output_path = sys.argv[2]

sn_path = pathlib.PurePath(input_path).parents[0]
sns = []
for item in os.listdir(sn_path):
    sn = item.split('.')[0]
    sns.append(sn)

current_sample = os.path.basename(input_path).split('.')[0]
pool = [current_sample]
pool += sample(sns, 99)

with open(output_path, 'w') as f:
    for item in pool:
        f.write('{0}\n'.format(item))