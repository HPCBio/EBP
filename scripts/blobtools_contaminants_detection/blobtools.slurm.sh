#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=300g
#SBATCH -N 1       #nodes
#SBATCH -n 12      #cpus-per-task
#SBATCH -p normal  #queue
#SBATCH -A ebp     #account
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@illinois.edu
#SBATCH -J EBxx_blobtools #jobname
#SBATCH -D /working/dir

### edit lines 11-13 above and variables below

### next lines are included to help troubleshoot the run
### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -x
set -euo pipefail
IFS=$'\n\t'
echo `hostname`

echo `date`
echo ################################ Declare variables and sanity check
### Declare variables
### These variables need editing
### GENOME is the genome after the sealer step with full path
### PREFIX is the prefix of output files
### OUTPUTDIR is the folder where results will be written to
### R1 and R2 are tellseq reads formatted for 10x
### Uniref should be diamind-indexed and its taxlist should also be available

PREFIX=EBxx
SCRATCH=/scratch/$SLURM_JOB_NAME
GENOME=/home/groups/earthbiogenome/results/some/path/Sealer/EBxx_Sealer_scaffold.fa
R1=/home/groups/earthbiogenome/data/TellSeq_like_Supernova_reads/EBxx_nobc_4tenx_S1_L001_R1_001.fastq.gz
R2=/home/groups/earthbiogenome/data/TellSeq_like_Supernova_reads/EBxx_nobc_4tenx_S1_L001_R2_001.fastq.gz
OUTPUTDIR=/home/groups/earthbiogenome/results/date_EBxx_blobtools/
blobdir=/home/groups/earthbiogenome/src/blobtools
TAXLIST=$blobdir/uniref100.taxlist
DIAMONDDB=$blobdir/uniref100.dmnd.gz
DBNAME=$(basename $DIAMONDDB)
GENOMENAME=$(basename $GENOME)
R1NAME=$(basename $R1)
R2NAME=$(basename $R2)
DAA=${PREFIX}_matches.daa
M8=${PREFIX}_matches.m8
TAGC=${PREFIX}_matches.daa.tagc
out_prefix=${PREFIX}_TellSeq_10x_format_trimmed
blob_prefix=${PREFIX}_blobplot

if [ ! -d "$blobdir" ]
then
	die "$blobdir path not found"
fi
if [ ! -s "$GENOME" ]
then
	die "$GENOME file not found"
fi
if [ ! -s "$R1" ]
then
	die "$R1 file not found"
fi
if [ ! -s "$R2" ]
then
	die "$R2 file not found"
fi
if [ ! -s "$DIAMONDDB" ]
then
	die "$DIAMONDDB file not found"
fi
if [ ! -s "$TAXLIST" ]
then
	die "$TAXLIST file not found"
fi

if [ ! -d "$SCRATCH" ]
then
	mkdir -p $SCRATCH
fi


echo `date`
echo ################################ Step 1: stage in files to scratch

cd  $SCRATCH
cp  $GENOME $SCRATCH
cp  $R1 $SCRATCH
cp  $R2 $SCRATCH
cp  $DIAMONDDB $SCRATCH 
cp  $TAXLIST $SCRATCH
cp  $SCRIPT $SCRATCH
gunzip uniref100.dmnd.gz

echo `date`
echo ################################ Step 2: run Diamond

module load DIAMOND/2.0.9-IGB-gcc-8.2.0
module load Perl/5.28.1-IGB-gcc-8.2.0



diamond blastx -p $SLURM_NPROCS -t $SCRATCH -d $DBNAME -q $GENOMENAME --evalue 1e-25 --long-reads --outfmt 100 -o $DAA 
exitcode=$?
if [ $exitcode -ne 0 ]
then
echo diamond blastx cmd failed
exit $exitcode
fi
if [ ! -s "$DAA" ]
then
	die "$DAA file not created"
fi

diamond view -a $DAA -f 6 -o $M8 #-a is view input file, -f 6 is tab-delimited format, -o is view output file

perl $blobdir/daa_to_tagc.pl $TAXLIST $DAA

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo daa_to_tagc.pl cmd failed
exit $exitcode
fi
if [ ! -s "$TAGC" ]
then
	die "$TAGC file not created"
fi


cp ${PREFIX}_matches* $OUTPUTDIR

echo `date`
echo ################################ Step 3: run bwa alignment w TellSeq reads

module purge
module load asset/1.0.3-IGB-gcc-8.2.0
module load minimap/2.18-IGB-gcc-8.2.0
module load BWA/0.7.17-IGB-gcc-8.2.0
module load SAMtools/1.11-IGB-gcc-8.2.0

bwa index $GENOMENAME

# TellSeq reads: align TellSeq reads to assembly fasta file with bwa; output is bam file
# 10x is part of asset and is used o reformat reads
# Usage: 10x [options...] <in1.fq/in1.fq.gz> <in2.fq/in2.fq.gz>

10x -l 18 -m 19 -c -p $out_prefix $R1NAME $R2NAME

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo 10x failed
exit $exitcode
fi

bwa mem -t $SLURM_NPROCS $GENOMENAME \
 ${out_prefix}_1.fq.gz \
 ${out_prefix}_2.fq.gz | samtools view -b - > ${out_prefix}_unsorted.bam

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo bwa mem cmd failed
exit $exitcode
fi

samtools sort -@ $SLURM_NPROCS -o ${out_prefix}_sorted.bam ${out_prefix}_unsorted.bam

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo samtools sort cmd failed
exit $exitcode
fi

samtools index ${out_prefix}_sorted.bam

cp ${out_prefix}* $OUTPUTDIR

echo `date`
echo ################################ Step 4: run blobtools

module purge
module load blobtools/1.1.1-IGB-gcc-4.9.4-Python-3.6.1

blobtools create -i $GENOMENAME -b ${out_prefix}_sorted.bam -t $TAGC -o $blob_prefix

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo blobtools create cmd failed
exit $exitcode
fi

blobtools plot -i ${blob_prefix}.blobDB.json --colours $blobdir/colors_phylum.txt -p 20
blobtools plot -i ${blob_prefix}.blobDB.json -r superkingdom --colours $blobdir/colors_Superkingdom.txt -p 20

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo blobtools plot cmd failed
exit $exitcode
fi

blobtools view -i ${blob_prefix}.blobDB.json -o ${blob_prefix}.blobDB.json_phylum 
blobtools view -i ${blob_prefix}.blobDB.json -o ${blob_prefix}.blobDB.json_superkingdom -r superkingdom

exitcode=$?
if [ $exitcode -ne 0 ]
then
echo blobtools view cmd failed
exit $exitcode
fi

cp ${blob_prefix}* $OUTPUTDIR

echo `date`
echo ################################ done. reset and exit

cd $outputdir
#rm -rf $SCRATCH
