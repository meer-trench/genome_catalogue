#!/bin/bash

#SBATCH --job-name=NAME
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

#module load miniconda3/4.10.3
#conda init bash
#conda activate MEER


