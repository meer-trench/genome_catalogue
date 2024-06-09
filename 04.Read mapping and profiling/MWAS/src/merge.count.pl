use strict;
my($REF,$DIR,$LV,$OUT)=@ARGV;

open R,"<$REF" or die "Cannot open $REF!\n";

my(%REF,%BINS,%SEQS, @samples);

while(<R>){
  chomp;
  my @F = split;
  $F[0] =~ /(.+\.\d+)\.(k\d+_\d+)/;
  my($bin,$kseq) = ($1,$2);
  $SEQS{$F[0]}{BIN} = $bin;
  $SEQS{$F[0]}{SEQ} = $kseq;
  $SEQS{$F[0]}{LEN} = $F[1];
  $BINS{$bin}{LEN} += $F[1];
}
# supply UNMAPPED
$SEQS{UNMAPPED}{BIN} = "UNMAPPED";
$SEQS{UNMAPPED}{SEQ} = "UNMAPPED";
$SEQS{UNMAPPED}{LEN} = 1;
$BINS{UNMAPPED}{LEN} = 1;

close R;

my @files = `ls $DIR/*.$LV.count`;
while(@files){
  my $file = shift @files;
  $file =~ /\/(.+).$LV.count/;
  my $sam = $1;
  push @samples, $sam;
  open I,"<$file";
  while(<I>){
    chomp;
    my @F = split;
    next unless exists $SEQS{$F[0]};
    my $bin = $SEQS{$F[0]}{BIN};
    my $len = $SEQS{$F[0]}{LEN};
    $SEQS{$F[0]}{$sam}{COUNT} = $F[1];
    #$SEQS{$F[0]}{$sam}{NORMC} = $F[1]*150 / $REF{$F[0]}{LEN};

    $BINS{$bin}{$sam}{COUNT} += $F[1];
    #$BINS{$bin}{$sam}{NORMC} += $SEQS{$F[0]}{$sam}{NORMC};
  }
  close I;
}

# Output bins profile
open B,">$OUT.bin.raw.count" or die "$!";
#open B,">$OUT.bin.norm.count" or die "$!";
print B join("\t","",@samples)."\n";
foreach my $b (sort keys %BINS){
  print B "$b";
  for(my $i=0; $i<@samples; $i++){
    my $count = ($BINS{$b}{$samples[$i]}{COUNT})?$BINS{$b}{$samples[$i]}{COUNT}:0;
    print B "\t$count";
  }
  print B "\n";
}
close B;

# Output all seqs profile
open S,">$OUT.seq.raw.count" or die "$!";
#open B,">$OUT.bin.norm.count" or die "$!";
print S join("\t","",@samples)."\n";
foreach my $s (sort keys %SEQS){
  print S "$s";
  for(my $i=0; $i<@samples; $i++){
    my $count = ($SEQS{$s}{$samples[$i]}{COUNT})?$SEQS{$s}{$samples[$i]}{COUNT}:0;
    print S "\t$count";
  }
  print S "\n";
}
close S;

exit;
