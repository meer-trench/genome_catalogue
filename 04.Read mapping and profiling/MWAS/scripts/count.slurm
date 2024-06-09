#!/bin/bash

#SBATCH --job-name=count
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/array_%A_%a.out
#SBATCH --error=tmp/array_%A_%a.err
#SBATCH --array=1-263%300

sh tmp/count.$SLURM_ARRAY_TASK_ID.sh
