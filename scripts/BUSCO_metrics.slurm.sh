#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J BUSCO_metrics #jobname
#SBATCH --mem=10g #memory requested
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -A ebp #account
#SBATCH -p normal #queue
#SBATCH --mail-user=your_email@illinois.edu #email address
#SBATCH --mail-type=ALL #get all job notifications
#SBATCH -D /home/groups/earthbiogenome/results/path/to/your/working/directory/

# ----------------Load Modules--------------------

module load Perl/5.28.1-IGB-gcc-8.2.0

module load BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2

# ----------------Commands------------------------


# modify your email in parameters above
# modify your working directory in parameters above

# change "hifasm_prefix" to correct file names


########## generate assembly metrics to go into shared Excel spreadsheet ###########

perl -I /home/groups/earthbiogenome/src/assemblathon2-analysis/ /home/groups/earthbiogenome/src/assemblathon2-analysis/assemblathon_stats.pl -csv hifiasm_prefix.p_ctg.fasta 

################## run BUSCO ##########################

export AUGUSTUS_CONFIG_PATH=/home/groups/earthbiogenome/src/config/

export BUSCO_CONFIG_FILE=/home/groups/earthbiogenome/src/config.ini

busco --offline -c 12 -i hifiasm_prefix.p_ctg.fasta -o BUSCO_hifiasm_prefix \
-l /home/apps/software/BUSCO/4.1.4-IGB-gcc-8.2.0-Python-3.7.2/data/lineages/insecta_odb10 -m genome -f 