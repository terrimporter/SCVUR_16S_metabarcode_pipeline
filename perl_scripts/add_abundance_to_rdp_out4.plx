#!/path/to/miniconda3/envs/myenv/bin/perl
# Teresita M. Porter, August, 21, 2019
# Script to add read abundance from ESV.table to rdp.out
# Prints to STDOUT so redirect to a file in snakemake
# USAGE perl add_abundance_to_rdp_out4.plx ESV.table rdp.out > outfile

use strict;
use warnings;

#declare var
my $i=0;
my $global_otu;
my $line;
my $assignment;
my $otu;
my $j=0;
my $sample;
my $abund;
my $outfile = "rdp.csv";
my $idline;
my $seqID;

#declare array
my @table;
my @rdp;
my @line;
my @headers;
my @otu;
my @idline;
my @sample;

#declare hash
my %assignment; #key = global OTU, value = rdp taxonomic assignment
my %table; #hash of hashes, key1=otu, key2=sample, value=readnumber

open (IN, "<", $ARGV[0]) || die "Error cannot open infile1: $!\n";
@table = <IN>;
close IN;

open (IN2, "<", $ARGV[1]) || die "Error cannot open infile2: $!\n";
@rdp = <IN2>;
close IN2;

#parse through global rdp taxonomic assignment, do not unmap the global OTUs here
while ($rdp[$i]) {
        $line = $rdp[$i];
        chomp $line;

        @line = split(/\t/,$line);
        $idline = shift @line;
        @idline = split(/;/, $idline);
        $global_otu = shift @idline;
        $assignment = join(",", @line); #csv output
        $assignment{$global_otu} = $assignment;
        $i++;
        @line=();

}
$i=0;

#grab list of headers and otus and build hash of hashes
while ($table[$i]) {
        $line = $table[$i];
        chomp $line;

        if ($i==0) { #header row
                @headers = split(' ',$line);#try splitting on whitespace instead of tabs
                shift @headers;
                shift @headers;
        }
        else {
                @line = split(' ',$line);#split on whitespace
                $global_otu = shift @line;
                $global_otu =~ s/^\s+//g; #strip any leading whitespace
                $global_otu =~ s/\s+$//g; #strip any trailing whitespace (do NOT use chomp)

                if (scalar @headers != scalar @line) {
                        print "Check FAILED\n";
                }

                foreach $abund (@line) {
                        $sample = $headers[$j];
                        $table{$global_otu}{$sample} = $abund;
                        $j++;
                }

                $j=0;
                $global_otu=();
                #$w++;
        }
        @line=();
        $i++;
}
$i=0;

#loop through hash of hashes, if OTU abund >= 3 keep and append taxonomic assignment, print out new assignment report
#open (OUT, ">>", $outfile) || die "Error cannot open outfile : $!\n";
print STDOUT "GlobalESV,SampleName,ESVsize,Strand,Root,RootRank,rBP,SuperKingdom,SuperKingdomRank,skBP,Kingdom,KingdomRank,kBP,Phylum,PhylumRank,pBP,Class,ClassRank,cBP,Order,OrderRank,oBP,Family,FamilyRank,fBP,Genus,GenusRank,gBP,Species,SpeciesRank,sBP\n";

#add read abundance to rdp out as a new column, unmap here!
while (($global_otu, $assignment) = each %assignment) {

        if (exists $table{$global_otu}) {

                while (($sample,$abund)  = each %{ $table{$global_otu}}) {
                        $abund = $table{$global_otu}{$sample};

                        if ($abund >= 3) {##### removes global OTUs with zero, one, or two mapped reads only #####

                                if ($assignment{$global_otu}) {
                                        $assignment = $assignment{$global_otu};
                                        print STDOUT "$global_otu,$sample,$abund,$assignment\n";
                                }
                                else {
                                        print "Assignment not available from RDP file\n";
                                }
                        }
                        else {
                                $j++;
                                next;
                        }
                }
        }
        else {
                print "Cannot find global_otu $global_otu in table\n";
        }
}

#close OUT;
