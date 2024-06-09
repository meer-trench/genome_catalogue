
#!/bin/bash

#SBATCH --job-name=jgi
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/jgi_%A_%a.out
#SBATCH --error=tmp/log/jgi_%A_%a.err
#SBATCH --array=8,132,138,142%4

sh tmp/jgi.$SLURM_ARRAY_TASK_ID.sh
