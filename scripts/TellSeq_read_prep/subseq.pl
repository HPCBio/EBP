#!/usr/bin/env perl
use Getopt::Long;
use strict;


######################################
# USAGE
######################################
my $USAGE =<<USAGE;

Script Name: subseq.pl 

Purposes: subseq.pl removes and reorders sequencing reads of R2 according to read order of R1
          
USAGE: perl subseq.pl -r1 inR1.fq -r2 inR2.fq -out ouR2.fq


USAGE

######################################
# parse input and sanity check
######################################

my $help;
my %R2=();
my $read1; my $read2; my $outfile;
my $commandLine1 = $0 . " " . (join " ", @ARGV);
my $n = @ARGV;

if (($commandLine1 =~ /-h/) || $n < 1){
   print "$USAGE\n";
   exit;
}

GetOptions ('r1=s'  => \$read1,
            'r2=s'  => \$read2,
            'out=s' => \$outfile,
            'h=s'   => \$help
);
            
if ( ($commandLine1 !~ /\s+-r1\s+/) || ($commandLine1 !~ /\s+-r2\s+/) || ($commandLine1 !~ /\s+-out\s+/) ) {
    print "Error in cmd: -r1 and -r2 -out are required\n\n";
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

open (OUT, ">", $outfile) or die "Can't open output file: $outfile\n";

######################################
# hash read2
######################################
my $totreads=0;
open inFile, $read2 or die "Can't open read2 file: $read2\n";
while(my $L1 = <inFile>){ 
	my $L2= <inFile>; my $L3= <inFile>; my $L4= <inFile>;
	my ($part1, $part2) = split(/ /,$L1);
	$R2{$part1}="$L1$L2$L3$L4";
 $totreads++;
}
close(inFile);
print "$read2 total reads: $totreads  \n";

######################################
# read1 compare and generate output
######################################
$totreads=0;
open inFile, $read1 or die "Can't open read1 file: $read1\n";
while(my $L1 = <inFile>){ 
	my $L2= <inFile>; my $L3= <inFile>; my $L4= <inFile>;
	my ($part1, $part2) = split(/ /,$L1);
        if ( defined $R2{$part1} ) {
		        print OUT "$R2{$part1}";
            $totreads++;       
        } else {
            print "ERROR: $part1 read not found in $read2. exiting now\n";
            exit 1;
        }
}
close(inFile);
print "$outfile total reads: $totreads  \n";

close(OUT);
%R2=();

      


