#!/path/to/miniconda3/envs/myenv/bin/perl
# By: Teresita M. Porter
# Date: August 20, 2019
# Script to add sample name from filenames to the fasta header
# Prints to STDOUT, in snakemake, redirect to outputfile
# USAGE perl rename_fasta_gzip.plx file.fasta.gz

use strict;
use warnings;

#declare vary
my $filename;
my $base;
my $base2;
my $outfile;
my $line;
my $i=0;
my $newline;
my $temp;
my $flag=0;

#declare array
my @in;
my @filename;

if ($ARGV[0] =~ /gz$/) {
	open (IN, "gunzip -c $ARGV[0] |") || die "Cannot open gzipped infile: $!\n";
	@in=<IN>;
	close IN;
	$flag=1;
}
else {
	open (IN, "<", $ARGV[0]) || die "Cannot open infile: $!\n";
	@in=<IN>;
	close IN;
}

$filename = $ARGV[0];
@filename = split(/\./,$filename); # split on dots in the filename
$base = $filename[0];
$base2 = (split '/', $base)[-1]; # split on forward slash to keep sample and exclude directory path, keep last element


while ($in[$i]) {
	$line = $in[$i];
	chomp $line;

	if ($line =~ /^>/) {
		$line =~ s/^>//g; 
		$newline = $base2.";".$line;
		print STDOUT ">".$newline."\n";
	}
	else {
		print STDOUT "$line\n";
	}
	$newline=();
	$i++;
}
$i=0;
