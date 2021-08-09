#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J purge_dups #jobname
#SBATCH --mem=50g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=youremail@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /path/to/your/working/directory/

# ----------------Load Modules--------------------
# For documentation, see https://github.com/dfguan/purge_dups

module load purge_dups/1.2.5-IGB-gcc-8.2.0-Python-3.7.2

# ----------------Commands------------------------

# change your email address above if using
# change your working directory above in -D, if using

# ---------------

# 1. create a text file in working directory called "pbfofn" listing the complete path(s) of the PacBio fastq file, one per line
# /home/groups/earthbiogenome/data/PacBio_HiFi_reads/HiFi.fastq

# ---------------

# 2. generate a config file for the purge_dups run, make sure you have correct path to reference fasta file
# usage: pd_config.py [-h] [-s SRF] [-l LOCD] [-n FN] [--version] ref pbfofn

# run this step interactively:

# srun -A ebp --pty bash

# module load purge_dups/1.2.5-IGB-gcc-8.2.0-Python-3.7.2

# ln -s /path/to/hifiasm.fasta .
# pd_config.py -l ./tmp_purge_hifiasm -n prefix_purge_hifiasm_config.json draft_assembly.fa pbfofn

# exit your interactive job

# ---------------

# 3. open .json config file and edit to skip BUSCO and KCM
# change busco and kcm blocks FROM "skip": 0  TO "skip": 1

# ---------------

# 4. run purge_dups w/config file
# usage: run_purge_dups.py [-h] [-p PLTFM] [-w WAIT] [-r RETRIES] [--version] config bin_dir spid

run_purge_dups.py -p bash prefix_purge_hifiasm_config.json /home/apps/software/purge_dups/1.2.5-IGB-gcc-8.2.0-Python-3.7.2/bin prefix_hifiasm_purge

# 5. Generate a plot of the cutoffs. Change the "reference_prefix" and" prefix_hist" to appropriate names.

hist_plot.py -c ./reference_prefix/coverage/cutoffs ./reference_prefix/coverage/PB.stat prefix_hist_plot.png