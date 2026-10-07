#!/usr/bin/env perl -w
use strict;

my $top = 10;

my (%matrix_75, %matrix_57, %rg, %rg_p7, %rg_p5, %all_p7, %all_p5, %all_rg, %cross, %cross2, %sum, %top, %expected_rg);

&help unless $ARGV[0];

open READ, $ARGV[0] or die "no input file was provided\n";
while (my $line = <READ>) {
    chomp $line;
    next unless $line =~ /^\d/;
    my ($number, $p7seq, $p7, $p5seq, $p5, $rg) = split /\t/, $line;
    next unless $p7 =~ /[^\-]/; #Removes sequences with unknown p7 index
    next unless $p5 =~ /[^\-]/; #Removes sequences with unknown p5 index
    next if $p7 =~ /phix[acgt]/i; #Removes sequences matching multi-base PhiX index (MPI-specific)
    next if $p5 =~ /phix[acgt]/i; #Removes sequences matching multi-base PhiX index (MPI-specific)
    $matrix_75{$p7}{$p5} = $number;
    $matrix_57{$p5}{$p7} = $number;
    if ($rg eq "unknown" || $rg eq "unexpected" || $rg eq "" || $rg eq "-") {
	$all_p5{$p5}++;
	$all_p7{$p7}++;
    } else {
	$expected_rg{$rg} = $number;
	$rg{$p7}{$p5} = $rg;
	$all_p5{$p5}++;
	$all_p7{$p7}++;
	$rg_p5{$p5}++;
	$rg_p7{$p7}++;
	$all_rg{$rg}{"p7"} = $p7;
	$all_rg{$rg}{"p5"} = $p5;
    }
}

foreach my $rg (keys %all_rg) {
    my $rg_p7 = $all_rg{$rg}{"p7"};
    my $rg_p5 = $all_rg{$rg}{"p5"};
    my $rg_count = $matrix_75{$rg_p7}{$rg_p5};
    foreach my $other_p5 (keys %{$matrix_75{$rg_p7}}) {
	my ($corner1, $corner2, $min_corner, $cont_name);
	next if $other_p5 eq $rg_p5;
	foreach my $other_p7 (keys %{$matrix_57{$other_p5}}) {
	    next if $other_p7 eq $rg_p7;
	    next unless exists $matrix_75{$other_p7}{$rg_p5};
	    $corner1 = $matrix_75{$rg_p7}{$other_p5};
	    $corner2 = $matrix_75{$other_p7}{$rg_p5};
	    $min_corner = $corner1;
	    $min_corner = $corner2 if $corner2 < $corner1;
	    my $cont = (($min_corner/$rg_count)**2 * $rg_count);
	    my $cont_name = $rg{$other_p7}{$other_p5} || "$other_p7/$other_p5";
	    $cross{"$cont_name into $rg"} = sprintf ("%0.3f", $cont);
	    $cross2{"$cont_name into $rg"} = $rg_count;
	    $sum{$rg} += $cont if $cont >= 0.5;
	    if (exists $top{$rg}{"num"}) {
		my $current = $top{$rg}{"num"};
		if ($cont > $current) {
		    $top{$rg}{"num"} = $cont if $cont > $current;
		    $top{$rg}{"source"} = $cont_name;
		}
	    }
	}
    }
}

my $counter = 0;
foreach my $element (sort {$cross{$b} <=> $cross{$a}} keys %cross) {
    $counter++;
    my $value = sprintf ("%.2f", $cross{$element});
    my $total = $cross2{$element};
    my $percent = sprintf ("%.5f", $value / $total * 100);
    print "#$element\t$value reads (of a total of ", $total, "), or $percent%\n";
    if ($counter >= $top) {
	last if $value < 0.5;
    }
}

print "\n#\#RG\tcross_cont_readsum\tcross_cont_percent\n";
foreach my $rg (sort keys %expected_rg) {
    print "$rg\t";
    if (exists $sum{$rg}) {
	my $total = $expected_rg{$rg};
	my $sum = sprintf ("%0.1f", $sum{$rg});
	print "$sum\t";
	my $percent = sprintf("%0.1f", $sum / $total * 1000000);
	print "$percent\n";
    } else {
	print "0\t0\n";
    }
}

print "#\n";

sub help {
    print "

This script requires a demultiplexing report to estimate possible cross-contamination between samples based on index swapping during library amplification or sequencing. 

The demultiplexing report has to be provided in the following format:
#number   p7seq     p7index   p5seq     p5index   readgroup
833843    CGATAGTT  298       TATGACGA  338       L4
693049    GTAGAATC  127       AGCATCGA  167       L1
583948    GGCCTGAC  9         CAACCTTA  49        L3
                         ...
2324      AGGTCGTC  -         CAACCTTA  49            <-- p7 index sequence unrecognized ('-'), readgroup unassigned
1839      ACTTAAGC  19        GAGGAGTC  372           <-- unexpected indexcombination, readgroup unassigned

Fields:
number    number of sequenced clusters, sorted by occurence of p7 and p5 index sequences
p7seq     observed p7 index sequence
p7index   Name/Number of index ('-' or empty if index is unrecognized)
p5seq     observed p5 index sequence
p5index   Name/Number of index ('-' or empty if index is unrecognized)
readgroup Library/sample name ('-','unknown','unexpected' or empty if index combination does not match an expected readgroup)

[usage]
./cross-contamination.pl demultiplexing-report.txt

[output]
STDOUT
1. List of all cross-contamination events with an estimated contribution of >0.5 reads (but at least the top 10 events, even if estimated to below 0.5 reads)
2. For each of the expected readgroups: sum of contaminating reads (only >0.5 per event), contribution of cross contamination in percent

";
exit;
}
