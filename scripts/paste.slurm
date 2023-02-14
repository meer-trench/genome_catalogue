#!/bin/bash

#SBATCH --job-name=paste
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/paste_%A.out
#SBATCH --error=tmp/log/paste_%A.err

sh tmp/paste.sh
