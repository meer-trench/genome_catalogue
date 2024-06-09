#!/bin/bash

#SBATCH --job-name=bowti2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --ntasks-per-node=64
#SBATCH --mail-type=end
#SBATCH --mail-user=fangchao@genomics.cn
#SBATCH --output=tmp/log/%j.out
#SBATCH --error=tmp/log/%j.err

kraken2 --db $K2DB_PF --threads 64 --report taxon/all_passed_bins.kreport2 drep_all/all_passed_bins.rename.fasta > taxon/all_passed_bins.kraken2 && echo "done"
