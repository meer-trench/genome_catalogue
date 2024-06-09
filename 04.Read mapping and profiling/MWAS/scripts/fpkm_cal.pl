#/usr/bin/env perl
use strict;

our ($insert_size,$STAT_FILE,$JGI_FILE) = @ARGV;
our (%STAT,%JGI);

open S,"<$STAT_FILE" or die $!;
while(<S>){
    chomp;
    $_ =~ /^([\S\s]+):(\s|\t)*(\d+)/;
    $STAT{$1} = $3;
}
close S;

my $N = $STAT{'raw total sequences'};
open J,"<$JGI_FILE" or die "cannot open $JGI_FILE ".$!;
$_=<J>;chomp; my @HEADS = split (/\t/,$_);
print $_."\tFPKM\n";
while(<J>){
    chomp;
    my @F = split (/\t/,$_);
    my $fpkm = $F[3] * 10**9 / ( $insert_size * $N );
    printf("%s\t%.4e\n",$_,$fpkm);
}
close J;

exit;