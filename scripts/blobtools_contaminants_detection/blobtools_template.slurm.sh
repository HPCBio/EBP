#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J blobtools  #jobname
#SBATCH --mem=200g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p hpcbio #queue
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/path/to/working/directory/

# ----------------Load Modules--------------------

module load DIAMOND/2.0.6-IGB-gcc-8.2.0
module load Perl/5.28.1-IGB-gcc-8.2.0

# ----------------Commands------------------------

# add uniref tax info to matches.daa file and format for blobtools
# creates matches.daa.tagc

perl daa_to_tagc.pl uniref100.taxlist matches.daa

module purge

# ----------------Load Modules--------------------

module load blobtools/1.1.1-IGB-gcc-4.9.4-Python-3.6.1

# ----------------Commands------------------------

# create the blobtools.json file

blobtools create -i supernova_fasta_file.fa -b prefix_sorted.bam -t matches.daa.tagc -o prefix_blobplot

# Creates plots

# Phylum level
blobtools plot -i prefix_blobplot.blobDB.json --colours colors_phylum.txt -p 20
# Superkindgdom level
blobtools plot -i prefix_blobplot.blobDB.json -r superkingdom --colours colors_Superkingdom.txt -p 20

# Creates tables

# Phylum level
blobtools view -i prefix_blobplot.blobDB.json -o prefix_phylum 
# Superkingdom level
blobtools view -i prefix_blobplot.blobDB.json -o prefix_superkingdom -r superkingdom
