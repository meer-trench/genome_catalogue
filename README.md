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




# derep_95_path
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
### Mapping fq to contig catalog
The script include following steps:  
1. Mapping fq to contig catalog  
2. sort sam and save as bam  
3. calculate mapping depth by using jgi method  

ref: https://github.com/RasmussenLab/vamb
```bash
perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1.fq/R2.fq/;
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/bowtie2.$array.sh");
  print FH "bowtie2 -p 30 -x drep_all/drep_all_95.rename -S bowtie2/$base.d95.sam -1 $F[0] -2 $p2 && \\\n"; 
  print FH "samtools view -@ 30 -b bowtie2/$base.d95.sam| samtools sort -@ 10 -o bowtie2/$base.d95.bam && \\\n";
  print FH "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth depth/$base.d95.raw.jgi bowtie2/$base.d95.bam \n";
  close FH;' tmp/tmp.R1.lst

sbatch scripts/bowtie2.slurm

```

#### hotfix
```bash

$ for i in `find depth/*.d95.raw.jgi -size -1000`;do b=`basename -s .d95.raw.jgi $i`;grep -n $b tmp/tmp.R1.lst;done|sed 's/:.*$//'|sort -n|tr '\n'
1,15,20,123,125,251,260,261,272,303,312,325,330,341,342,349,351,360,363,367,375,380,403,410,411,433,436,452,460,468,513,515,523,535,539,541,552,566,576,591,608,615,618,625,659,679,695,775,785,799,805,838,863,869,919,927,949,953,960,984,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,1036,1037,1038,1039,1040,1041,1042,1043,1044,1045,1046,1047,1048,1049,1050,1051,1052,1053,1054,1055,1056,1057,1058,1059,1060,1061,1062,1063,1064,1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,1115,1116,1117,1118,1119,1120,1145,1149,1150,1153,1167
#manually add 695 for FDZ071-RRR2-4.1


```

```bash

find bowtie2/*.d95.bam -size -1000|sed 's/.d95.bam//;s/bowtie2\///'|sort > empty.bam.lst
find depth/*.d95.raw.jgi -size -1000|sed 's/.d95.raw.jgi//;s/depth\///'|sort > empty.jgi.lst
comm -23 empty.jgi.lst empty.bam.lst > valid.bam.lst

for i in `cat valid.bam.lst`;do b=`grep -n $i tmp/tmp.R1.lst|sed 's/:.*$//'`;sed -i '1s/^/#/' tmp/bowtie2.$b.sh;done

cp scripts/bowtie2.slurm scripts/bowtie2.patch.slurm
# Add above arrays into slurm file
sbatch scripts/bowtie2.patch.slurm
```

### jgi to PFKM
After fnish all jgi calculation. try to generate a PFKM profile.
Use bam-readcount to generate reads count information from bam files: https://github.com/genome/bam-readcount
```bash
perl -ane '($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1.fq/R2.fq/;
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++;
  open(FH, "> tmp/stat.$array.sh");
  print FH "samtools view -h bowtie2/$base.d95.bam|samtools stats|grep ^SN | cut -f 2- > bowtie2/$base.d95.stat \n";
  close FH;' tmp/tmp.R1.lst

sbatch scripts/stat.slurm
```

patch to unify sample naming：
```bash
perl -ane 'BEGIN{print "sample\tbase\told_base\n"};($base=`basename $F[0]`)=~s/_R1.fq.gz//;
  chomp($base);
  $old = $base;
  $base =~ s/_$//g;$base =~ s/_/-/g;
  $sample = $base;
  if(exists $HS{$old}){ $old .= ".1"; $base .= ".1";}; $HS{$old}++;
  if(exists $BS{$base}){ $base .= ".1";}; $BS{$base}++;
  print "$sample\t$base\t$old\n"' tmp/tmp.R1.lst > tmp/base.name.lst
sed -i 's/^FDZ064-BuG10-12-7/FDZ064-BuG10-12/' base.name.lst

perl -ane 'chomp; $old=$F[2]; $base=$F[1]; $sample=$F[0];
  `mv depth/$old.d95.MAG.fpkm depth/$base.d95.MAG.fpkm`;
  ' base.name.lst

perl -ane 'chomp; next if $F[0] eq "sample";
  $old=$F[2]; $base=$F[1]; $sample=$F[0]; $array++;
  open(FH, "> tmp/fpkm.$array.sh");
  print FH "perl scripts/fpkm_cal.pl 150 bowtie2/$old.d95.stat depth/$old.d95.raw.jgi > depth/$base.d95.fpkm\n";
  close FH;' tmp/base.name.lst

sbatch scripts/fpkm.slurm
```

```bash
perl -ane 'chomp; $old=$F[2]; $base=$F[1]; $sample=$F[0];
  $array++; $s++;$l=1; 
  open(FH, "< depth/$base.d95.fpkm");
  while(<FH>){chomp;@FS=split; $HS{V}{$FS[0]}{$base} = $FS[4]; 
    $HS{R}{$FS[0]}{C} ++ if $FS[4] > 0; $HS{R}{$FS[0]}{L} = $FS[1];
    $HS{C}{$base}{C} ++ if $FS[4] > 0;  $HS{C}{$base}{S} += $FS[4]; $l++;
  }; close FH; $M=$l if $l > $M; print STDERR "processed sample(#$s) \r"; 
  END{print STDERR "\nSummarizing ... \r"; open R,">profiles/MEER1225.fpkm.row.stat"; 
    open C,">profiles/MEER1225.fpkm.col.stat"; open V,">profiles/MEER1225.fpkm.profile";
    print R "OTU\tlength\tcounts\n"; print C "sample\tsum\tcounts\n";
    print V "OTU"; foreach $c (sort keys %{$HS{C}}){
      print V "\t$c"; print C "$c\t$HS{C}{$c}{S}\t$HS{C}{$c}{C}\n";
    }; print V "\n"; $l=0;
    foreach $r (sort keys %{$HS{R}}){
      print R "$r\t$HS{R}{$r}{L}\t$HS{R}{$r}{C}\n"; print V "$r"; $l++;
      foreach $c (sort keys %{$HS{C}}){
        print V "\t$HS{V}{$r}{$c}";
      }
      print V "\n";
      print STDERR "Summarizing $l / $M \r";
    };print STDERR "\nDONE\n"
  }' base.name.lst

```
This will generate `profiles/MEER1225.fpkm.profile` wiht FPKM of all 1225 samples. The FPKM was calculated from coverage by using `jgi_summarize_bam_contig_depths`. This is a method from Maxbat2.

$$
jgi_i = \frac{X_i * ins} {l_i}
$$

where $X_i$ is read counts mapped to contig $i$; $ins$ is insert size, $l_i$ is conig effective length (trimmed ends).

$$
X_i = \frac {jgi \cdot l_i}  {ins}
$$

FPKM: Substitute reads with fragments per kilobase of exon/gene/contig/genome per million reads mapped.

$$
FPKM = \frac {X_i} {(l_i / 10^3)(N/10^6)} = \frac {X_i} {l_iN} \cdot 10^9
$$

$$
=\frac {jgi_i\cdot l_i / ins} {l_iN} \cdot 10^9 = \frac {jgi_i} {N \cdot ins} \cdot 10^9
$$

Where $N$ is total number of reads sequenced.

‍Get MAG profile:

```bash
perl -ane 'chomp; $old=$F[2]; $base=$F[1]; $sample=$F[0];
  $array++; $HS{$base}++; $s++;$l=1; 
  open(FH, "> tmp/MAG95.$array.sh");
  print FH "perl scripts/MAG_cal.pl 150 profiles/MEER1225.fpkm.row.stat.gz drep_all/drep_all_95.Widb.csv depth/$base.d95.fpkm > depth/$base.d95.MAG.fpkm\n";
  close FH;' base.name.lst

sbatch scripts/fpkm.MAG95.slurm

#merge to be a single table
perl -ane 'BEGIN{%HS;$s=0;$l=0;$M=1;};
  ($base=`basename $F[0]`)=~s/_R1.fq.gz//;chomp($base); ($p2=$F[0])=~s/R1.fq/R2.fq/;
  if(exists $HS{$base}){$base .= ".1"}; $array++; $HS{$base}++; $s++;$l=1; 
  open(FH, "< depth/$base.d95.MAG.fpkm");<FH>;
  while(<FH>){chomp;@FS=split; $HS{V}{$FS[0]}{$base} = $FS[4]; 
    $HS{R}{$FS[0]}{C} ++ if $FS[4] > 0; $HS{R}{$FS[0]}{L} = $FS[1];
    $HS{C}{$base}{C} ++ if $FS[4] > 0;  $HS{C}{$base}{S} += $FS[4]; $l++;
  }; close FH; $M=$l if $l > $M; print STDERR "processed sample(#$s) \r"; 
  END{print STDERR "\nSummarizing ... \r"; open R,"|bgzip -c >profiles/MEER1225.d95.fpkm.MAG.row.stat.gz"; 
    open C,"|bgzip -c >profiles/MEER1225.d95.fpkm.MAG.col.stat.gz"; open V,"|bgzip -c >profiles/MEER1225.d95.fpkm.MAG.profile.gz";
    print R "OTU\tlength\tcounts\n"; print C "sample\tsum\tcounts\n";
    print V "OTU"; foreach $c (sort keys %{$HS{C}}){
      print V "\t$c"; print C "$c\t$HS{C}{$c}{S}\t$HS{C}{$c}{C}\n";
    }; print V "\n"; $l=0;
    foreach $r (sort keys %{$HS{R}}){
      print R "$r\t$HS{R}{$r}{L}\t$HS{R}{$r}{C}\n"; print V "$r"; $l++;
      foreach $c (sort keys %{$HS{C}}){
        print V "\t$HS{V}{$r}{$c}";
      }
      print V "\n";
      print STDERR "Summarizing $l / $M \r";
    };print STDERR "\nDONE\n"
  }' tmp/tmp.R1.lst

ls profiles/MEER1225.d95.fpkm.MAG.*|grep -E -v "dvc|gz"|xargs -n1 bgzip
dvc add profiles/MEER1225.d95.fpkm.MAG.*.gz
```

vamb:
```bash
vamb -p 60 --outdir VAMB --fasta drep_all/all_passed_bins.rename.fasta --jgi jgi/jgi.abundance.dat 
```

```bash

vamb --outdir VAMB10 --fasta drep_all/all_passed_bins.10.fasta --jgi jgi/jgi.abundance.10.dat 

```
### Adding fungi profile from Zewei
```bash
perl -e 'open I,"< profiles/fungi.jgi.depth.tsv";
$head = <I>; chomp($head); @heads = split(/\t/,$head);
for ($i=0;$i<$#heads;$i++){ if($heads[$i] =~ /sorted.bam$/){
  $heads[$i] =~ s/.sorted.bam//; 
  $heads[$i] =~ s/_//; 
  $HS{$i} = $heads[$i];
}};
while(<I>){
  chomp;
  @FS = split /\t/;
  $INDEX{$FS[0]}{len} = $FS[1];
  $INDEX{$FS[0]}{avg} = $FS[2];
  foreach $i(sort {$a<=>$b} keys %HS){
    $HS{$i}{$FS[0]}{dep} = $FS[3];
    $HS{$i}{$FS[0]}{var} = $FS[3];
  }
}
foreach $i(sort {$a<=>$b} keys %HS){
  open I,"> depth/$HS{$i}.fungi.jgi";
  print I "contigName\tcontigLen\ttotalAvgDepth\t$HS{$i}.fungi.bam\t$HS{$i}.fungi.bam-var\n";
  foreach $c(sort keys %INDEX){
    print I "$c\t$INDEX{$c}{len}\t$INDEX{$c}{avg}\t";
    print I "$HS{$i}{$c}{dep}\t$HS{$i}{$c}{var}\n";
  }
  close I;
}'

mkdir bowtie2/sample.stat
perl -ane 'chomp; $old=$F[2]; $base=$F[1]; $sample=$F[0];
  open I,"<bowtie2/$old.d95.stat";
  while(<I>){chomp;@FS=split /:\s+/;$HS{$sample}{$FS[0]} +=$FS[1]};
  close I;
  END{
    foreach $s(keys %HS){
      open O,">bowtie2/sample.stat/$s.d95.stat";
      foreach $k (keys %{$HS{$s}}){
        print O "$k: $HS{$s}{$k}\n";
      }
      close O;
    }
  }
' tmp/base.name.lst

perl -ane 'chomp; next if $F[0] eq "sample"; $old=$F[2]; $base=$F[1]; $sample=$F[0]; $array++;
  `perl scripts/fpkm_cal.pl 150 bowtie2/sample.stat/$sample.d95.stat depth/$sample.fungi.jgi > depth/$base.fungi.fpkm`;
' tmp/base.name.lst

# get MAG profile

perl -ane 'chomp; next if $F[0] eq "sample"; $old=$F[2]; $base=$F[1]; $sample=$F[0];
  next if exists $HS{C}{$sample};
  $array++; $s++;$l=1;
  open(FH, "< depth/$sample.fungi.fpkm");
  while(<FH>){chomp;@FS=split; next if $FS[0] eq "contigName";
    $HS{V}{$FS[0]}{$sample} = $FS[5]; 
    $HS{R}{$FS[0]}{C} ++ if $FS[5] > 0; $HS{R}{$FS[0]}{L} = $FS[1];
    $HS{C}{$sample}{C} ++ if $FS[5] > 0;  $HS{C}{$sample}{S} += $FS[5]; $l++;
  }; close FH; $M=$l if $l > $M; print STDERR "processed sample(#$s) \r"; 
  END{print STDERR "\nSummarizing ... \r"; 
  open R,"|bgzip -c >profiles/MEER1194.fungi.fpkm.row.stat.gz"; 
  open C,"|bgzip -c >profiles/MEER1194.fungi.fpkm.col.stat.gz";
  open V,"|bgzip -c >profiles/MEER1194.fungi.fpkm.profile.gz";
    print R "OTU\tlength\tcounts\n"; print C "sample\tsum\tcounts\n";
    print V "OTU"; foreach $c (sort keys %{$HS{C}}){
      print V "\t$c"; print C "$c\t$HS{C}{$c}{S}\t$HS{C}{$c}{C}\n";
    }; print V "\n"; $l=0;
    foreach $r (sort keys %{$HS{R}}){
      print R "$r\t$HS{R}{$r}{L}\t$HS{R}{$r}{C}\n"; print V "$r"; $l++;
      foreach $c (sort keys %{$HS{C}}){
        print V "\t$HS{V}{$r}{$c}";
      }
      print V "\n";
      print STDERR "Summarizing $l / $M \r";
    }; print STDERR "\nDONE\n";
  }' tmp/base.name.lst

dvc add profiles/MEER1194.fungi.fpkm.{{col,row}.stat,profile}.gz

```
### Taxonomic classiy

```bash
perl -e 'my $MAP_FILE="drep_all/drep_all_95.Widb.csv";
open M, ($MAP_FILE =~ /.gz$/)?"bgzip -dc $MAP_FILE|":"<$MAP_FILE" or die $!;
my $MAP_header = <M>;
chomp($MAP_header);
my @MAP_heads = split(/,/,$MAP_header);
my @cluster_index = grep {$MAP_heads[$_] eq "cluster"}  0..$#MAP_heads;
while(<M>){chomp; my @F = split /,/;
    $F[0] =~ /^(\S+)\.bin\.(\S+)\.fa$/;
    print "$1.$2\t$F[$cluster_index[0]]\n";
}' > drep_all/drep_all_95.binmap.tsv
```

Try Zewei's binners
```bash
DIR=/dssg/home/acct-trench/trench/USER/songzewei/data_process/tst_binning_pathways


/dssg/home/acct-trench/trench/USER/songzewei/data_process/tst_binning_pathways/data/derep_pathways/data_tables/Wdb.csv

/dssg/home/acct-trench/trench/USER/songzewei/data_process/tst_binning_pathways/data/derep_pathways/dereplicated_genomes/FDZ064-BaGW12-14.bin.107.fa|

perl -e ' open I, "<drep_all/drep99_gtdbtk_classify/meer_drep99_all.summary.tsv";
$C=0; <I>; while(<I>){
  chomp; @F=split /\t/; @F0 = split(/.bin./,$F[0]); $pfx = "$F0[0].$F0[1]";
  $HS{$pfx} = $C; $C++;
}; close I; print STDERR "Done file1.\n";
open L,"<drep_all/drep_all_99.rename.lst";
while (<L>){
  chomp; @F= split /\./; $pfx="$F[0].$F[1]";
  if(exists $HS{$pfx}){
    print "$HS{$pfx}\t$_\n";
  }
};close L;' |bgzip > taxon/drep99_gtdbtk_classify.clusters.lst.gz

for i in {S,G,F,O,C,P};do 
perl scripts/binTaxonBench.pl $i taxon/drep99_gtdbtk_classify.clusters.lst.gz \
  taxon/all_passed_bins.kraken2.detail.gz \
  > taxon/drep_99.gtdbtk.$i.benchmark.tsv &
done

#d99
perl -e ' open I, "< drep_all/drep_all_99/data_tables/Wdb.csv";
  <I>; while(<I>){
  chomp; @F=split /,/; @F0 = split(/\./,$F[0]); $pfx = "$F0[0].$F0[2]";
  $HS{$pfx} = $F[1];
}; close I; print STDERR "Done file1.\n";
open L,"<drep_all/drep_all_99.rename.lst";
while (<L>){
  chomp; @F= split /\./; $pfx="$F[0].$F[1]";
  if(exists $HS{$pfx}){
    print "$HS{$pfx}\t$_\n";
  }
};close L;' |bgzip > taxon/drep99_cluster.lst.gz

for i in {S,G,F,O,C,P};do 
perl scripts/binTaxonBench.pl $i taxon/drep99_cluster.lst.gz \
  taxon/all_passed_bins.kraken2.detail.gz \
  > taxon/drep_99.cluster.$i.benchmark.tsv &
done
```


# test profile
```bash
perl scripts/call.prof.pl S jgi/jgi.abundance.dat taxon/all_passed_bins.kraken2.detail.gz > profiles/jgi.kk2.S.prof

```