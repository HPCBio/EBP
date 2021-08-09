#!/bin/bash

# ----------------SLURM Parameters----------------

#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email@illinois.edu
#SBATCH -J EBxx_SALSA #jobname
#SBATCH -D /path/to/working/dir

# ----------------Load Modules--------------------

module load BEDTools/2.28.0-IGB-gcc-8.2.0

# ----------------Commands------------------------

# Change email and working path directory above as needed.
# Change paths and prefixes below as needed.

# Move into arima "deduplicated_files" directory and convert bam to bed for input to SALSA.

cd /home/groups/earthbiogenome/results/assembly/path/arima_out/deduplicated_files/

bedtools bamtobed -i prefix_rep1.bam > prefix_rep1.bed

# Sort bedfile by read name:

sort -k 4 prefix_rep1.bed > prefix_rep1_name_sort.bed

module purge

# ----------------Load Modules--------------------

module load SAMtools/1.11-IGB-gcc-8.2.0

# ----------------Commands------------------------

# SAMtools index genome assembly:

cd /home/groups/earthbiogenome/results/assembly/

samtools faidx prefix_hifiasm_purged_scaffolded_assembly.fa

module purge

# ----------------Load Modules--------------------

module load SALSA/20191001-IGB-gcc-4.9.4-Python-2.7.13

# ----------------Commands------------------------

# Run SALSA

cd /home/groups/earthbiogenome/results/20210305_purged_scaffolded_EB19/salsa_out/

run_pipeline.py -a ../prefix_hifiasm_purged_scaffolded_assembly.fa \
 -l ../prefix_hifiasm_purged_scaffolded_assembly.fa.fai \
 -b ../arima_out/deduplicated_files/prefix_rep1_name_sort.bed \
 -e DNASE -o scaffolds -m yes -p yes
 
#usage: run_pipeline.py [-h] -a PATH TO ASSEMBLY FASTA FILE -l PATH TO ASSEMBLY FAI FILE -b PATH TO BED FILE [-o OUTPUT DIRECTORY NAME]
#                       [-c CUTOFF] [-g GFA] [-u UNITIGS] [-e ENZYME: for OMNI-C = DNASE]
#                       [-i ITER] [-x DUP] [-s EXP] [-m CLEAN/CORRECT ERRORS] [-p output the scaffolds sequence and agp file for each iteration]
