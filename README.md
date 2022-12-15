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
Pick following represented samples:
FDZ077-WuYW0-2
FDZ077-WuYW2-4
FDZ077-WuYW4-6
FDZ077-WuYW8-10
FDZ077-WuYW12-14
FDZ064-RRG0-2
FDZ064-RRG2-4
FDZ064-RRG4-6
FDZ064-RRG8-10
FDZ064-RRG12-14
FDZ050-RBk0-2
FDZ050-RBk2-4
FDZ050-RBk4-6
FDZ050-RBk8-10
FDZ050-RBk12-14
FDZ060-HbG0-2
FDZ060-HbG2-4
FDZ060-HbG4-6
FDZ060-HbY0-2
FDZ060-HbY2-4
FDZ060-HbY4-6
FDZ060-HbY8-10
FDZ060-HbY12-14
FDZ036-RY0-2
FDZ036-RY2-4
FDZ036-RY4-6
FDZ036-RY8-10
FDZ036-RY12-14
FDZ041-Be0-5
FDZ041-Be5-10
FDZ041-Be10-15
FDZ033-YBk0-2
FDZ033-YBk2-4
FDZ033-YBk4-6
FDZ033-YBk8-10
FDZ033-YBk12-14
FDZ033-GG0-2
FDZ033-GG2-4
FDZ033-GG4-6
FDZ033-GG8-10
FDZ033-GG12-14
FDZ036-Y0-2
FDZ036-Y2-4
FDZ036-Y4-6
FDZ036-Y8-10
FDZ036-Y12-14

prepare reference:
```bash
for i in `cat select.46.sample.name.lst|xargs -n1`;do find /dssg/home/acct-trench/share/DATA/RAW_Merge -name "$i*R1.fq.gz";done > tmp/tmp.R1.lst

for i in `cat tmp.R1.lst`;do b=$(echo $i| tr R1 R2);n=`basename $i|sed 's/_R1.fq.gz//' `;sed 's/NAME/$n/g;s/PATH1/$i/;s/PATH2/$b/' tmp/temp.slurm > tmp/$n.slurm;done

conda activate meer
# mkdir bowtie2

# build index
srun -N 1 -n 1 -p 64c512g bowtie2-build --threads 60 drep_all/drep_all_95.rename.fasta drep_all/drep_all_95.rename


```

mapping:
```bash
perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1/R2/;
  `sed "s/NAME/$base/" temp.slurm > tmp/$base\.slurm`; open(FH, ">> tmp/$base.slurm");
  print FH "bowtie2 -p 30 -x drep_all/drep_all_95.rename -S bowtie2/$base.d95.sam \\
  -1 $F[0] \\\n  -2 $p2\n"; close FH;' tmp/tmp.R1.lst

# run
for i in `ls tmp/*.slurm`;do sbatch sbatch $i; done
```

Count mapping reads:
```bash
for i in  `ls bowtie2/*.sam|xargs -n1 basename|sed 's/.d95.sam//'`; do
  srun -N 1 -n 1 -p 64c512g perl src/counts.pl bowtie2/$i.d95.sam bowtie2/$i.d95.count &
done
```

Sumamrize reference length:
```bash
module load samtools/1.13-gcc-11.2.0
samtools faidx drep_all/drep_all_95.rename.fasta

srun -N 1 -n 1 -p 64c512g perl src/merge.count.pl drep_all/drep_all_95.rename.fasta.fai bowtie2 d95 select33.d95.sum
```

Visualization can be found in `VISUAL.Rmd`.
