#/usr/bin/env perl
use strict

# Get taxon rank tree from kraken2 db and then 
# assigned them to the contigs
our($INSPECT_FILE,$SEQLIST_FILE) = @ARGV;

# Rank order: R D K P C O F G S
our %RANKS = (
    'S' => 90,
    'G' => 80,
    'F' => 70,
    'O' => 60,
    'C' => 50,
    'P' => 40,
    'K' => 30,
    'D' => 20,
    'R' => 10,    
);
foreach my $i (R,D,K,P,C,O,F,G,S) {
    foreach my $s ("",1,2,3,4,5,6,7,8,9){
        my $rank  = "$i$s";
        my $value = $RANKS{$i} + $s;
        $RANKS{$rank} = $value;
        $RANKS{$value} = $rank;
    }
}

my (%R_PARENTS,%TMP_PARENTS);
my $preRankValue = 999;
my $preRankID    = 0;

open INS, "<$INSPECT_FILE" or die $!;

while(<INS>){
    next if /^#/;
    chomp;
    my @Fs = split /\t/;
    $Fs[3] =~ /^(\S)(|\d)$/;
    my ($mainRank,$subRank) = ($1, $2);
    my $curRankValue = $RANKS{$mainRank} + $subRank;
    for(my $i=$curRankValue; $i<=99;$i++){
        delete $TMP_PARENTS{$RANKS{$i}};
    };
    $TMP_PARENTS{$Fs[3]} = $Fs[4];
    %{$R_PARENTS{$Fs[4]}} = %TMP_PARENTS;

}
close INS;

open SEQ,"bgzip -dc $SEQLIST_FILE|" or die $!;
while (<SEQ>){
    chomp;
    my @Fs = split /\t/;
    my ($S,$G,$F,$O,$C,$P,$K,$D,$R) = (
        $R_PARENTS{$Fs[2]}{S}?$R_PARENTS{$Fs[2]}{S}:0,
        $R_PARENTS{$Fs[2]}{G}?$R_PARENTS{$Fs[2]}{G}:0,
        $R_PARENTS{$Fs[2]}{F}?$R_PARENTS{$Fs[2]}{F}:0,
        $R_PARENTS{$Fs[2]}{O}?$R_PARENTS{$Fs[2]}{O}:0,
        $R_PARENTS{$Fs[2]}{C}?$R_PARENTS{$Fs[2]}{C}:0,
        $R_PARENTS{$Fs[2]}{P}?$R_PARENTS{$Fs[2]}{P}:0,
        $R_PARENTS{$Fs[2]}{K}?$R_PARENTS{$Fs[2]}{K}:0,
        $R_PARENTS{$Fs[2]}{D}?$R_PARENTS{$Fs[2]}{D}:0,
        $R_PARENTS{$Fs[2]}{R}?$R_PARENTS{$Fs[2]}{R}:0
    );
    printf ("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$_,$S,$G,$F,$O,$C,$P,$K,$D,$R);
}

exit;