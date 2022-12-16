#

# Quick glance
## merge kraken2 profiles
```bash
perl src/get.kk2.prf.pl P
perl src/get.kk2.prf.pl G
```
Visualization can be found in `VISUAL.Rmd`.

# Generate MAGs profile
Fowllowing codes were initially tested in SJTU HPC:

Using SRUN: https://docs.hpc.sjtu.edu.cn/job/slurm.html#sbatch
```bash
# Interactive method:
salloc -N 1 -n 16 -p 64c512g
ssh casxxx
```

```bash
ln -s ../../songzewei/data_process ../

mkdir drep_all
indir="../data_process/drep_data/drep_all_95/dereplicated_genomes"
for i in `ls $indir`;do
  awk -v f=$i 'FNR==1{sub(".bin","",f);sub("fa","",f)};/^>/{sub(">",">"f,$0)}{print}' $indir/$i
done > drep_all/drep_all_95.rename.fasta


```


prepare reference:
```bash
# build index
conda activate meer
srun -N 1 -n 1 -p 64c512g bowtie2-build --threads 60 drep_all/drep_all_95.rename.fasta drep_all/drep_all_95.rename
```

find samples from selected dives:
```bash
for i in {FDZ081,FDZ077,FDZ064,FDZ050,FDZ060,FDZ036,FDZ041,FDZ033};do find /dssg/home/acct-trench/share/DATA/RAW_Merge -name "$i*R1.fq.gz";done > tmp/tmp.R1.lst
```
found 284 samples, but 277 unique.

# supp.
for i in {FDZ061,FDZ062,FDZ063};do find /dssg/home/acct-trench/share/DATA/RAW_Merge -name "$i*R1.fq.gz";done > tmp/tmp.R1.lst
## add 87 samples

mapping:
```bash
mkdir bowtie2

perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1.fq/R2.fq/;
  next if exists $HS{$base}; $array++; $HS{$base}++;
  open(FH, "> tmp/bowtie2.$array.sh");
  print FH "bowtie2 -p 30 -x drep_all/drep_all_95.rename -S bowtie2/$base.d95.sam \\
  -1 $F[0] \\\n  -2 $p2\n"; close FH;' tmp/tmp.R1.lst

# submit working batch
conda activate meer
sbatch scripts/bowtie2.slurm
```

Count mapping reads:
```bash
j=1
for i in  `ls bowtie2/*.sam|xargs -n1 basename|sed 's/.d95.sam//'`; do
  echo perl src/counts.pl bowtie2/$i.d95.sam bowtie2/$i.d95.count > tmp/count.$j.sh;
  j=$(($j + 1));
done
# alternatively
j=1
for i in  `find bowtie2 -name "*.sam" -atime -0.25|xargs -n1 basename|sed 's/.d95.sam//'`; do
  echo perl src/counts.pl bowtie2/$i.d95.sam bowtie2/$i.d95.count > tmp/count.$j.sh;
  j=$(($j + 1));
done

#
sbatch scripts/count.slurm

```

Sumamrize reference length:
```bash
module load samtools/1.13-gcc-11.2.0
samtools faidx drep_all/drep_all_95.rename.fasta

salloc -N 1 -n 40 -p 64c512g
ssh nodexxx
nohup perl src/merge.count.pl drep_all/drep_all_95.rename.fasta.fai bowtie2 d95 profiles/select11dives.d95.sum &
```

Visualization can be found in `VISUAL.Rmd`.
