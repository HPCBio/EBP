#!/usr/bin/perl -w

use strict;
use warnings;
# use this script to convert an Illumina fastq file to fasta format.
# fastq format has 4 lines: @sequence header, sequence, +sequence header, quality score
# read in file line by line; if line begins with '@', store as $line1, print to outfile
# read next line and print it.
# skip two lines
# repeat

my $filein = $ARGV[ 0 ];
open IN, $filein or die "Cannot open infile\n";

my $out_file = $filein . '.fasta';
open OUT, ">$out_file" or die "Cannot open outfile\n";

my $header;
my $sequence;

while ( $header = <IN> ){ 
	chomp $header;
	$header=~/^@(\S+)/;
	print OUT ">$1\n";
	$sequence = <IN>;
	print OUT $sequence;
	<IN>; 
	<IN>;
}

close IN;
close OUT;
exit;