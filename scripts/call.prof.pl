#/usr/bin/env perl
use strict;


my ($RANK, $RAW_PROF_FILE, $ANNO_FILE) = @ARGV;

our (%ANNOMAP,%PROF);

our %RPOS = (
    "S" => 4, "G" => 5, "F" => 6, "O" => 7, "C" => 8,
    "P" => 9, "K" => 10, "D" => 11, "R" => 12
);

open A, "bgzip -dc $ANNO_FILE|" or die $!;
while(<A>){
    chomp;
    my @F= split /\t/;
    $ANNOMAP{$F[1]} = $F[$RPOS{$RANK}];
}
close A;
print STDERR "Reading raw profile ... ";

open P, "<$RAW_PROF_FILE" or die $!;
my $HEAD = <P>;
my @HEADS = split (/\t/,$HEAD);
my @SAMPLES = @HEADS;
shift @SAMPLES; shift @SAMPLES; shift @SAMPLES;
my $counts =0;
while(<P>){
    chomp;
    my @F = split /\t/;
    my $len = $F[1];
    my $taxonID = $ANNOMAP{$F[0]};
    for(my $i=3;$i<=$#HEADS;$i++){
        my $sample = $HEADS[$i];
        $PROF{$taxonID}{$sample}{weighted_depth} += $len * $F[$i];
        $PROF{$taxonID}{$sample}{weights} += $len;
    }
    $counts++;
    if($counts % 10000 == 0){
        print STDERR "Reading raw profile ... lines: $counts \r";
    }
    
}
print STDERR "\nReading raw profile ... Done.\n Summarizing ... ";

print join("\t","taxonID",@SAMPLES)."\n";
foreach my $taxonID (sort {$a <=> $b} keys %PROF){
    print "taxonID";
    for(my $i=0;$i<=$#SAMPLES;$i++){
        my $sample = $SAMPLES[$i];
        my $count = $PROF{$taxonID}{$sample}{weighted_depth} / $PROF{$taxonID}{$sample}{weights};
        print "\t$count";
    }
    print "\n";
}
close P;

print STDERR "Done!\n";

exit;

