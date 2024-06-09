#/usr/bin/env perl
use strict;

our ($RANK,$CLUSTER_FILE,$SEQTAXON_FILE,$DREP_LIST) = @ARGV;
our (%ID2CLUST,%CLUSTINFO,%GENOMEINFO,%DREP);

if($DREP_LIST){
    print STDERR "Reading $DREP_LIST ...";
    open DREP, "< $DREP_LIST" or die $!;
    while(<DREP>){
        chomp;
        $DREP{$_} = 1;
    }
    close DREP;
    print STDERR " Done.\n";
}

print STDERR "Reading $CLUSTER_FILE ...";

open CLUSTER, "bgzip -dc $CLUSTER_FILE|" or die $!;
while(<CLUSTER>){
    chomp;
    my @F = split /\t/;
    $ID2CLUST{$F[1]} = $F[0];
}
close CLUSTER;

print STDERR " Done.\nReading $SEQTAXON_FILE ...";

open TAXON, "bgzip -dc $SEQTAXON_FILE|" or die $!;
my %RANKPOS = (
    S =>  4, G =>  5, F =>  6, O => 7, C => 8, P => 9,
    K => 10, D => 11, R => 12
);
while(<TAXON>){
    chomp;
    my @F = split /\t/;
    next if $F[0] eq "U";
    if($DREP_LIST){
        next unless $DREP{$F[1]} == 1;
    }
    my $cluster = ($ID2CLUST{$F[1]})?$ID2CLUST{$F[1]}:0;
    my $taxon = $F[$RANKPOS{$RANK}];
    $CLUSTINFO{$cluster}{L}{$F[3]} += $F[3];
    $CLUSTINFO{$cluster}{LS}       += $F[3];
    $CLUSTINFO{$cluster}{$taxon}   += $F[3];

    $GENOMEINFO{$taxon}{$cluster}  += $F[3];
    $GENOMEINFO{$taxon}{SUM} += $F[3];

}
close TAXON;
print STDERR " Done.\nSummarizing ...";

# Calculate
foreach my $c (sort {$a<=>$b} keys %CLUSTINFO){
    my $members = keys %{$CLUSTINFO{$c}{L}};
    my $sumlength = $CLUSTINFO{$c}{LS};

    my ($taxon, $TP, $FP, $FN, $purity, $completeness) = &TBFPFN($c);

    #
    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\n",$c,$members,$sumlength,$taxon, $TP, $FP, $FN, $purity, $completeness);
    
}

print STDERR "Done. All Done!\n";
exit;

#########################################
# TP: True positive
sub TBFPFN{
    my $c = shift;
    my ($g,$TP,$FP,$FN) = (0, $CLUSTINFO{$c}{0},$CLUSTINFO{$c}{LS} - $CLUSTINFO{$c}{0},$GENOMEINFO{0}{SUM} - $GENOMEINFO{0}{$c});
    foreach my $ig (sort { $CLUSTINFO{$c}{$b} <=> $CLUSTINFO{$c}{$a} } keys %{$CLUSTINFO{$c}}) {
        if($ig !=0 && $g == 0){
            $g = $ig;
            $TP = $CLUSTINFO{$c}{$g};
            $FN = $GENOMEINFO{$g}{SUM} - $GENOMEINFO{$g}{$c};
            $FP = $CLUSTINFO{$c}{LS} - $CLUSTINFO{$c}{$g};
            last;
        }
    }

    my $purity = ($TP==0)?0: $TP / ($TP + $FP);
    my $completeness = ($TP==0)?0: $TP / ($TP + $FN);
    return($g, $TP, $FP, $FN, $purity, $completeness);
}