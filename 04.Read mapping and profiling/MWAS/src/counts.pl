use strict;
my($i,$o) = @ARGV;

open I,"<$i" or die "Cannot open $i!\n";
open O,">$o" or die "Cannot open $o!\n";

my (%REF);

while(<I>){
  next if $_=~/^@/;
  chomp;
  my @F = split;
  if($F[2] eq "*" || $F[4] eq "*"){
    $REF{"UNMAPPED"} ++;
  }else{
    $REF{$F[2]} ++;
  }
}
close I;

foreach my $r (sort keys %REF){
  print O "$r\t$REF{$r}\n";
}
close O;
exit;
