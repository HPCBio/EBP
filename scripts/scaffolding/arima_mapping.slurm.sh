#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH -J arima_mapping #jobname
#SBATCH -D /home/groups/earthbiogenome/results/assembly/path/arima_out/

# ----------------Load Modules--------------------

module load SAMtools/1.11-IGB-gcc-8.2.0 

module load BWA/0.7.17-IGB-gcc-8.2.0

module load picard/2.10.1-Java-1.8.0_152

# ----------------Commands------------------------

# Change email and working path directory above as needed.
# Change sample information below as needed.


##############################################
# ARIMA GENOMICS MAPPING PIPELINE 02/08/2019 #
##############################################

#Below find the commands used to map HiC data.

#Replace the variables at the top with the correct paths for the locations of files/programs on your system.

#This bash script will map one paired end HiC dataset (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with your volume of samples and system.

##########################################
# Commands #
##########################################

#SRA='basename_of_fastq_files'
#LABEL='overall_exp_name'
#BWA='/path/to/bwa'
#SAMTOOLS='/path/to/samtools'
#IN_DIR='/path/to/gzipped/fastq/files'
#REF='/path/to/reference_sequences.fa'
#FAIDX='$REF.fai'
#PREFIX='bwa_index_name'
#RAW_DIR='/path/to/write/out/bams'
#FILT_DIR='/path/to/write/out/filtered/bams'
#FILTER='/path/to/filter_five_end.pl'
#COMBINER='/path/to/two_read_bam_combiner.pl'
#STATS='/path/to/get_stats.pl'
#PICARD='/path/to/picard.jar'
#TMP_DIR='/path/to/write/out/temporary/files'
#PAIR_DIR='/path/to/write/out/paired/bams'
#REP_DIR='/path/to/where/you/want/deduplicated/files'
#REP_LABEL=$LABEL\_rep1
#MERGE_DIR='/path/to/final/merged/alignments/from/any/biological/replicates'
#MAPQ_FILTER=10
#CPU=12

SRA='EB19_CAGATC_L003'
LABEL='EB19_hifiasm_purged_scaff'
BWA='/home/apps/software/BWA/0.7.17-IGB-gcc-8.2.0/bin/bwa'
SAMTOOLS='/home/apps/software/SAMtools/1.11-IGB-gcc-8.2.0/bin/samtools'
IN_DIR='/home/groups/earthbiogenome/data/HiC_reads'
REF='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/EB19_hifiasm_purged_scaffolded_prettyNames.fa'
FAIDX='$REF.fai'
PREFIX='EB19_hifiasm_purged_scaffolded_prettyNames'
RAW_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/raw_bams'
FILT_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/filtered_bams'
FILTER='/home/groups/earthbiogenome/src/filter_five_end.pl'
COMBINER='/home/groups/earthbiogenome/src/two_read_bam_combiner.pl'
STATS='/home/groups/earthbiogenome/src/get_stats.pl'
PICARD='/home/apps/software/picard/2.10.1-Java-1.8.0_152/picard.jar'
TMP_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/tmp_files'
PAIR_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/paired_bams'
REP_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/deduplicated_files'
REP_LABEL=$LABEL\_rep1
MERGE_DIR='/home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/arima_out/replicates'
MAPQ_FILTER=10
CPU=12


echo "### Step 0: Check output directories exist & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
$BWA index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
$BWA mem -t $CPU $PREFIX $IN_DIR/$SRA\_R1_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
$BWA mem -t $CPU $PREFIX $IN_DIR/$SRA\_R2_001.fastq.gz | $SAMTOOLS view -@ $CPU -Sb - > $RAW_DIR/$SRA\_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_1.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
$SAMTOOLS view -h $RAW_DIR/$SRA\_2.bam | perl $FILTER | $SAMTOOLS view -Sb - > $FILT_DIR/$SRA\_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/$SRA\_1.bam $FILT_DIR/$SRA\_2.bam $SAMTOOLS $MAPQ_FILTER | $SAMTOOLS view -bS -t $FAIDX - | $SAMTOOLS sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $PICARD AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=ILLUMINA PU=none

###############################################################################################################################################################
###                                           How to Accommodate Technical Replicates                                                                       ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                                                    ###
### Technical replicates (eg. one library split across multiple lanes) should be merged before running the MarkDuplicates command.                          ###
### If this step is run, the names and locations of input files to subsequent steps will need to be modified in order for subsequent steps to run correctly.###
### The code below is an example of how to merge technical replicates.                                                                                      ###
###############################################################################################################################################################
#	REP_NUM=X #number of the technical replicate set e.g. 1
#	REP_LABEL=$LABEL\_rep$REP_NUM
#	INPUTS_TECH_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as technical replicates
#   example bash array - INPUTS_TECH_REPS=('INPUT=A.L1.bam' 'INPUT=A.L2.bam' 'INPUT=A.L3.bam')
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_TECH_REPS OUTPUT=$TMP_DIR/$REP_LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $PICARD MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

$SAMTOOLS index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"

#########################################################################################################################################
###                                       How to Accommodate Biological Replicates                                                    ###
### This pipeline is currently built for processing a single sample with one read1 and read2 fastq file.                              ###
### Biological replicates (eg. multiple libraries made from the same sample) should be merged before proceeding with subsequent steps.###
### The code below is an example of how to merge biological replicates.                                                               ###
#########################################################################################################################################
#
#	INPUTS_BIOLOGICAL_REPS=('bash' 'array' 'of' 'bams' 'from' 'replicates') #BAM files you want combined as biological replicates
#   example bash array - INPUTS_BIOLOGICAL_REPS=('INPUT=A_rep1.bam' 'INPUT=A_rep2.bam' 'INPUT=A_rep3.bam')
#
#	java -Xmx8G -Djava.io.tmpdir=temp/ -jar $PICARD MergeSamFiles $INPUTS_BIOLOGICAL_REPS OUTPUT=$MERGE_DIR/$LABEL.bam USE_THREADING=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT
#
#	$SAMTOOLS index $MERGE_DIR/$LABEL.bam

# perl $STATS $MERGE_DIR/$LABEL.bam > $MERGE_DIR/$LABEL.bam.stats

# echo "Finished Mapping Pipeline through merging Biological Replicates"

