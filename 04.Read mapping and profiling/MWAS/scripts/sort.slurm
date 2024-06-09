#!/bin/bash

#SBATCH --job-name=patch
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/sort_%A_%a.out
#SBATCH --error=tmp/log/sort_%A_%a.err
#SBATCH --array=8,132,138,142%4

sh tmp/sort.$SLURM_ARRAY_TASK_ID.sh
