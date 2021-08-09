#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J SALSA_Juicebox_hic_assembly_files #jobname
#SBATCH --mem=100g
#SBATCH -N 1 #nodes
#SBATCH -n 8 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/assembly/path/salsa_out/scaffolds/

# ----------------Load Modules--------------------

module load SAMtools/1.11-IGB-gcc-8.2.0

# ----------------Commands------------------------

# Change email and working path directory above as needed.
# Change paths and prefixes below as needed.

# Index final SALSA fasta file (creates a tab-delimited file, first two columns are scaffold name and length):

samtools faidx scaffolds_FINAL.fasta

#Reformat .fai file to keep only first two columns.

cut -f 1,2 scaffolds_FINAL.fasta.fai > chromosome_sizes.tsv

module purge

# ----------------Load Modules--------------------

module load SALSA/20191001-IGB-gcc-4.9.4-Python-2.7.13

# ----------------Commands------------------------

#Convert alignments to text using bed file, agp file, and scaffold length file.

alignments2txt.py -b alignment_iteration_1.bed -a scaffolds_FINAL.agp -l scaffold_length_iteration_1 > alignments.txt

# '-b','--bed',help='Original bed file of alignments'
# '-a','--agp',help='AGP file for the final scaffolds'
# '-l','--length',help='Length of input unitigs'

#Sort the alignments.

awk '{if ($2 > $6) {print $1"\t"$6"\t"$7"\t"$8"\t"$5"\t"$2"\t"$3"\t"$4} else {print}}' alignments.txt | sort -k2,2d -k6,6d -T ./ --parallel=8 | awk 'NF' > alignments_sorted.txt

#Run Juicer "pre" to generate the .hic file for Juicebox. The .hic file describes the contacts relationships of the reads to reference fasta file.
#I ran this step several times to find the minimum amount of heap space; otherwise java fails.

java -Xmx20g -jar /home/apps/software/Juicer/1.6.0-IGB-gcc-8.2.0/juicer_tools.jar pre alignments_sorted.txt salsa_scaffolds_FINAL.hic chromosome_sizes.tsv

# Usage: <infile>: Text file with paired contacts, <outfile>: Name of outfile, should end with .hic, <genomeID>: path of the chromosome_sizes.tsv file

#====================================================

#Create the ".assembly" file for viewing in Juicebox, specific format that describes relationships of scaffolds:

awk -f /home/groups/earthbiogenome/src/generate-sorted-assembly-file-from-fasta.awk scaffolds_FINAL.fasta > scaffolds_FINAL.assembly