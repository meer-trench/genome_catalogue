#!/bin/bash

#SBATCH --job-name=cut
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/cut_%A_%a.out
#SBATCH --error=tmp/log/cut_%A_%a.err
#SBATCH --array=1-1225%100

sh tmp/cut45.$SLURM_ARRAY_TASK_ID.sh
