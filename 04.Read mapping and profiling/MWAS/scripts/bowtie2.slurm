#!/bin/bash

#SBATCH --job-name=bowtie2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/array_%A_%a.out
#SBATCH --error=tmp/array_%A_%a.err
#SBATCH --array=1-1225%256

sh tmp/bowtie2.$SLURM_ARRAY_TASK_ID.sh
