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
    $BINS{$F[$cluster_index[0]]} = "$1.$2";
}
close M;

open S, ($STAT_FILE =~ /.gz$/)?"bgzip -dc $STAT_FILE|":"<$STAT_FILE" or die $!;
while(<S>){
    chomp;
    my @F = split /\t/;
    $F[0] =~ /^(\S+)\.(\S+)\.(\S+)$/ or next;
    $MAP{"$1.$2"}{$3} = $F[1];
    $MAP{"$1.$2"}{"_sumLen"} += $F[1];
    $MAP{"$1.$2"}{"_contigs"} ++;
}
close S;

my $N = $STAT{'raw total sequences'};
open P,"<$PROFILE" or die $!;
$_=<P>;chomp; my @HEADS = split (/\t/,$_);
while(<P>){
    chomp;
    my @F = split (/\t/,$_);
    $F[0] =~ /^(\S+)\.(\S+)\.(\S+)$/ or next;
    $MAP{"$1.$2"}{"_sumFPKMxLEN"} += $F[1] * $F[4];
    $MAP{"$1.$2"}{"_sumLen_counted"} += $F[1];
    $MAP{"$1.$2"}{"_contigs_counted"} ++;
}
close P;

# summary
$HEADS[2] = "contigs";
$HEADS[3] = "theo_contigs";
print join("\t",@HEADS)."\n";
foreach my $bin (sort keys %BINS){
    my $b = $BINS{$bin};
    my $fpkm = 0;
    if($MAP{$b}{"_sumLen_counted"}>0){
        $fpkm = $MAP{$b}{"_sumFPKMxLEN"} / $MAP{$b}{"_sumLen_counted"};
    }
    if($MAP{$b}{"_contigs_counted"} ne $MAP{$b}{"_contigs"}){
        print STDERR "$bin\($b\) mapped $MAP{$b}{_contigs} contigs, ";
        print STDERR "but in practice $MAP{$b}{_contigs_counted} counted.\n";
    }
    if($MAP{$b}{"_sumLen_counted"} ne $MAP{$b}{"_sumLen"}){
        print STDERR "$bin\($b\) summed $MAP{$b}{_sumLen} bases, ";
        print STDERR "but in practice $MAP{$b}{_sumLen_counted} summped.\n";
    }
    printf("%s\t%d\t%d\t%d\t%.4e\n",$bin,$MAP{$b}{"_sumLen_counted"},
    $MAP{$b}{"_contigs_counted"},$MAP{$b}{"_contigs"},$fpkm);
}
exit;