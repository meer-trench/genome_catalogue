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




# derep_all_path
```bash
cd  /dssg/home/acct-trench/share/DATA/CLEAN/

indir="/dssg/home/acct-trench/trench/USER/songzewei/data_process/drep_data/all_passed_bins"
for i in `ls $indir`;do
  awk -v f=$i 'FNR==1{sub(".bin","",f);sub("fa","",f)};/^>/{sub(">",">"f,$0)}{print}' $indir/$i
done > all_passed_bins.rename.fasta
```

```bash
cd meer_sediments_metagenomics/mwas

cp /dssg/home/acct-trench/share/DATA/CLEAN/all_passed_bins.rename.fasta drep_all/
# build index
conda activate meer
srun -N 1 -n 40 -p 64c512g
salloc -N1 -n 24 -p huge
ssh xxx
bowtie2-build --threads 40 drep_all/all_passed_bins.rename.fasta drep_all/all_passed_bins.rename

```

```bash
rsync -atvP trench-1@pilogin:/lustre/home/acct-trench/trench-1/USER/fangchao/repo/meer_sediments_metagenomics/mwas/drep_all/all_passed_bins.rename.{1,2,3,4,rev.1,rev.2}.bt2l drep_all/
find /dssg/home/acct-trench/share/DATA/RAW_Merge -name "*R1.fq.gz" > tmp/tmp.R1.lst


bowtie2 -p 30 -x drep_all/all_passed_bins.rename -S bowtie2/FDZ071-RRR10-12.all.sam \
-1/dssg/home/acct-trench/share/DATA/RAW_Merge/sample627_0914/FDZ071-RRR10-12_R1.fq.gz -2 /dssg/home/acct-trench/share/DATA/RAW_Merge/sample627_0914/FDZ071-RRR10-12_R1.fq.gz
```

```bash
perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1.fq/R2.fq/;
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/bowtie2.$array.sh");
  print FH "bowtie2 -p 30 -x drep_all/all_passed_bins.rename -S bowtie2/$base.all.sam \\
  -1 $F[0] \\\n  -2 $p2\n"; close FH;' tmp/tmp.R1.lst


```
https://github.com/RasmussenLab/vamb

### Sort sam to bam
```bash
perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base);
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/sort.$array.sh");
  print FH "samtools view -@ 10 -b bowtie2/$base.all.sam| samtools sort -@ 10 -o bowtie2/$base.all.bam \n"; close FH;' tmp/tmp.R1.lst

sbatch scripts/sort.slurm
```

Calculate bins depth:
```bash
jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input} 2>{log}

perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base);
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/jgi.$array.sh");
  print FH "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth depth/$base.all.raw.jgi bowtie2/$base.all.bam \n"; close FH;' tmp/tmp.R1.lst

sbatch scripts/jgi.slurm
```

cut column 1-3, 4-5:
```bash
cut -f1-3 depth/FDZ031YE0-5.all.raw.jgi > jgi/jgi.column1to3

perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base);
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/cut45.$array.sh");
  print FH "cut -f1-3 --complement depth/$base.all.raw.jgi > depth/$base.all.4to5.jgi \n"; close FH;' tmp/tmp.R1.lst
sbatch scripts/cut.slurm

#paste jgi/jgi.column1to3 depth/*.all.4to5.jgi  > jgi/jgi.abundance.dat
# too many open files

ls depth/*.all.4to5.jgi|split -d -l 250 - tmp/split.jgi.

for i in {0..4};do
  cat tmp/split.jgi.0$i|xargs paste > tmp/paste.$i &
done
wait && paste jgi/jgi.column1to3 tmp/paste.{0..4} > jgi/jgi.abundance.dat

sbatch scripts/paste.slurm

```

vamb:
```bash
vamb -p 60 --outdir VAMB --fasta drep_all/all_passed_bins.rename.fasta --jgi jgi/jgi.abundance.dat 
```

```bash

vamb --outdir VAMB10 --fasta drep_all/all_passed_bins.10.fasta --jgi jgi/jgi.abundance.10.dat 

```