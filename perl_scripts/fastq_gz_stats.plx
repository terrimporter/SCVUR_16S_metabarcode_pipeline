#!/path/to/miniconda3/envs/myenv/bin/perl
# Aug. 13/19 edit to work with snakemake
# Teresita M. Porter
# Script to grab seq stats from fastq.gz files
# Write to STDOUT so redirect to an outfile with snakemake
# USAGE perl fastq_stats.plx infile.fastq.gz > outfile

use strict;
use warnings;
use Statistics::Lite qw(:all);

#declare var
my $num;
my $i=0;
my $line;
my $length;
my $num2;
my $count;
my $min;
my $max;
my $mean;
my $median;
my $mode;

#declare array
my @allseqs;
my @line;
my @lengths;

@allseqs = `zcat $ARGV[0] | awk '(NR%4==2)'`;
$num = scalar (@allseqs);

while ($allseqs[$i]) {
	$line = $allseqs[$i];
	chomp $line;

	@line = split(//,$line);
	$length = scalar @line;
	push @lengths, $length;
	$i++;
}
$i=0;

$num2 = scalar (@lengths);

if ($num != $num2) {
	print "Possible error\n";
}

$count = count (@lengths);
$min = min (@lengths);
$max = max (@lengths);
$mean = mean (@lengths);
$median = median (@lengths);
$mode = mode (@lengths);

#open (OUT, ">>", $ARGV[1]) || die "Error cannot open outfile:$!\n";

print STDOUT $ARGV[0]."\t".$count."\t".$min."\t".$max."\t".$mean."\t".$median."\t".$mode."\n";

#close OUT;

$i=0;
