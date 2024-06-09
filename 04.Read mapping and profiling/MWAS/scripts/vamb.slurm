#!/bin/bash

#SBATCH --job-name=vamb
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 60
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/vamb_%A.out
#SBATCH --error=tmp/log/vamb_%A.err

vamb -p 60 --outdir VAMB --fasta drep_all/all_passed_bins.rename.fasta --jgi jgi/jgi.abundance.dat
