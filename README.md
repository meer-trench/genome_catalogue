#


## merge kraken2 profiles

```perl
use strict;

our ($HASH,@SAMPLES);

my $R = $1;
my @FILES = split("\n", `ls \$MEER_SJTU/trench_MEER_v1_database_krakan2_diversity/*.report.txt`);

open P ">merge.$R.pct";
open C ">merge.$R.sct";
open X ">merge.$R.idx";

foreach my $f (@FILES){
  $f =~ /_diversity\/(.*).report.txt/;
  my $sid = $1;
  push @SAMPLES, $sid;
  print P "\t$sid";
  print C "\t$sid";

  open I,"<$f";
  while(<I>){
    chomp;
    my @F = split;
    next unless $F[3] =~ /^(U|$R)$/;
    $F[5] =~ s/\s+//g;
    $HASH{$F[4]}{rank} = $F[3];

    $HASH{$F[4]}{$sid}{pct} = $F[0];
    $HASH{$F[4]}{$sid}{sct} = $F[1];
    $HASH{$F[4]}{$sid}{ect} = $F[2];
    $HASH{$F[4]}{$sid}{annot} = $F[5];
  }
  close I;
}
print P "\n";
print C "\n";

foreach my $t (sort {$a<=>$b} keys %HASH){
  print X "$t\t$HASH{$t}{rank}\n";
  print P "$t";
  print C "$t";
  
  for (my $i=0; $i<=$#SAMPLES; $i ++){
    my $pct = ($HASH{$t}{$SAMPLES[$i]}{pct})?$HASH{$t}{$SAMPLES[$i]}{pct}:0;
    my $sct = ($HASH{$t}{$SAMPLES[$i]}{sct})?$HASH{$t}{$SAMPLES[$i]}{sct}:0;
    print P "\t$pct";
    print C "\t$sct";
  }
  print P "\n";
  print C "\n";
}
close P;
close C;
close X;

exit;

```

```bash
perl get.kk2.prf.pl P
perl get.kk2.prf.pl G
```
