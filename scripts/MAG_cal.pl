#/usr/bin/env perl
use strict;

our ($insert_size,$STAT_FILE,$MAP_FILE,$PROFILE) = @ARGV;
our (%STAT,%MAP,%PROF,%BINS);

open M, ($MAP_FILE =~ /.gz$/)?"bgzip -dc $MAP_FILE|":"<$MAP_FILE" or die $!;
my $MAP_header = <M>;
chomp($MAP_header);
my @MAP_heads = split(/,/,$MAP_header);
my @cluster_index = grep {$MAP_heads[$_] eq "cluster"}  0..$#MAP_heads;
while(<M>){
    chomp;
    my @F = split /,/;
    $F[0] =~ /^(\S+)\.bin\.(\S+)\.fa$/;
    $MAP{"$1.$2"}{"_cluster"} = $F[$cluster_index[0]];
    $MAP{$F[$cluster_index[0]]}{"$1.$2"} ++;
}
close M;

open S, ($STAT_FILE =~ /.gz$/)?"bgzip -dc $STAT_FILE|":"<$STAT_FILE" or die $!;
while(<S>){
    chomp;
    my @F = split /\t/;
    $F[0] =~ /^(\S+)\.(\S+)\.(\S+)$/ or next;
    my $cluster = $MAP{"$1.$2"}{"_cluster"};
    $STAT{$cluster}{"_sumLen"} += $F[1];
    $STAT{$cluster}{"_contigs"} ++;
}
close S;

open P,"<$PROFILE" or die $!;
$_=<P>;chomp; my @HEADS = split (/\t/,$_);
while(<P>){
    chomp;
    my @F = split (/\t/,$_);
    $F[0] =~ /^(\S+)\.(\S+)\.(\S+)$/ or next;
    my $cluster = $MAP{"$1.$2"}{"_cluster"};
    $BINS{$cluster}{"_sumFPKMxLEN"} += $F[1] * $F[4];
    $BINS{$cluster}{"_sumLen"} += $F[1];
    $BINS{$cluster}{"_contigs"} ++;
}
close P;

# summary
$HEADS[2] = "contigs";
$HEADS[3] = "theo_contigs";
print join("\t",@HEADS)."\n";
foreach my $b (sort keys %BINS){
    my $fpkm = 0;
    if($BINS{$b}{"_sumLen"}>0){
        $fpkm = $BINS{$b}{"_sumFPKMxLEN"} / $BINS{$b}{"_sumLen"};
    }
    if($BINS{$b}{"_contigs"} ne $STAT{$b}{"_contigs"}){
        print STDERR "$b mapped $STAT{$b}{_contigs} contigs, ";
        print STDERR "but in practice $BINS{$b}{_contigs} counted.\n";
    }
    if($BINS{$b}{"_sumLen"} ne $STAT{$b}{"_sumLen"}){
        print STDERR "$b summed $STAT{$b}{_sumLen} bases, ";
        print STDERR "but in practice $BINS{$b}{_sumLen} summped.\n";
    }
    printf("%s\t%d\t%d\t%d\t%.4e\n",$b,$BINS{$b}{"_sumLen"},
    $BINS{$b}{"_contigs"},$STAT{$b}{"_contigs"},$fpkm);
}
exit;