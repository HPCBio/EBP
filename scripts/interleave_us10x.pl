#!/usr/bin/env perl
use Getopt::Long;
use strict;


######################################
# USAGE
######################################
my $USAGE =<<USAGE;

Script Name: interleave_us10x.pl 

Purpose:
It reads tellseq reads and produces a single file of 10x-formatted interleaved reads that can be accepted by ARCS
It embeds barcodes from index file in each header of each read of R1 and R2
It interleaves R1 and R2. Finaly, it produces one file  

USAGE: perl interleave_us10x.pl -r1 read1.fq -r2 read2.fq -index index.fq -out interleaved.fq

All parameters are required

USAGE

######################################
# parse input and sanity check
######################################

my $help;
my %R2=(); my %Idx=();
my $read1; my $read2; my $outfile; 
my $index; my $read1_idx; my $read2_idx;

my $commandLine1 = $0 . " " . (join " ", @ARGV);
my $n = @ARGV;

if (($commandLine1 =~ /-h/) || $n < 1){
   print "$USAGE\n";
   exit;
}

GetOptions ('r1=s'  => \$read1,
            'r2=s'  => \$read2,
            'out=s' => \$outfile,
            'index=s' => \$index,            
            'h=s'   => \$help
);
            
if ( ($commandLine1 !~ /\s+-r1\s+/) || ($commandLine1 !~ /\s+-r2\s+/) ) {
    print "Error in cmd: -r1 and -r2 -out are required\n\n";
    exit 1;
}
if ( ($commandLine1 !~ /\s+-index\s+/) || ($commandLine1 !~ /\s+-out\s+/) ) {
    print "Error in cmd: -index -out are required\n\n";
    exit 1;
}

if (! -s $read1) {
    print "Error in cmd: -r1 $read1 file not found\n\n";
    exit 1;
}

if (! -s $read2) {
    print "Error in cmd: -r2 $read2 file not found\n\n";
    exit 1;
}

if (! -s $index) {
    print "Error in cmd: -index $index file not found\n\n";
    exit 1;
}

open (OUT, ">", $outfile) or die "Can't open output file: $outfile\n";

######################################
# hash index
######################################

open inFile, $index or die "Can't open index file: $index\n";
while(my $L1 = <inFile>){ 
	my $L2= <inFile>; my $L3= <inFile>; my $L4= <inFile>;
	my ($part1, $part2) = split(/ /,$L1);
	$Idx{$part1}="$L2";
}
close(inFile);

######################################
# hash read2
######################################

open inFile, $read2 or die "Can't open read2 file: $read2\n";
while(my $L1 = <inFile>){ 
	my $L2= <inFile>; my $L3= <inFile>; my $L4= <inFile>;
	my ($part1, $part2) = split(/ /,$L1);
	if ( defined $Idx{$part1} ) {
		$read1_idx = $Idx{$part1};
    $L1 = $part1; chomp $L1;
		$L1 = $L1." BX:Z:".$read1_idx;
		$R2{$part1}="$L1$L2$L3$L4";
	}
}
close(inFile);

######################################
# read1 compare and generate output
######################################

open inFile, $read1 or die "Can't open read1 file: $read1\n";
while(my $L1 = <inFile>){ 
	my $L2= <inFile>; my $L3= <inFile>; my $L4= <inFile>;
	my ($part1, $part2) = split(/ /,$L1);
	if ( defined $Idx{$part1} ) {
		$read1_idx = $Idx{$part1};
    $L1 = $part1; chomp $L1;
		$L1 = $L1." BX:Z:".$read1_idx;
        	if ( defined $R2{$part1} ) {
			 print OUT "$L1$L2$L3$L4";
			 my $read2_val = $R2{$part1};
			 print OUT "$read2_val";
		 } else { print "skipped read, could not find it in R2  $part1 \n"; }
	} else { print "skipped read, could not find it in Index  $part1 \n"; }
}
close(inFile);
close(OUT);
%R2=();
%Idx=();
      


