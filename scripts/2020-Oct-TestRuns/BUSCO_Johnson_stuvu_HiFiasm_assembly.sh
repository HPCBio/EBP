#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J 20201029_BUSCO_stuvu_Hifiasm
#SBATCH --mem=50g
#SBATCH -N 1 
#SBATCH -n 8 
#SBATCH -p normal
#SBATCH --mail-user=grendon@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D  /home/groups/hpcbio/projects/CBC/2020-Oct-HiCanu/results/2020-10-19-Hifiasm/
#SBATCH -o 20201029_BUSCO_stuvu_Hifiasm.ou
#SBATCH -e 20201029_BUSCO_stuvu_Hifiasm.er

# ----------------Load Modules--------------------

module load BUSCO/3.0.1-IGB-gcc-4.9.4-Python-2.7.13

export AUGUSTUS_CONFIG_PATH=/home/a-m/grendon/augustus/config

# ----------------Commands------------------------

cd /home/groups/hpcbio/projects/CBC/2020-Oct-HiCanu/results/2020-10-19-Hifiasm/

run_BUSCO.py -i stuvu_assembly.asm.p_ctg.fa -o BUSCO_stuvu_hifiasm \
-l /home/apps/software/BUSCO/3.0.1-IGB-gcc-4.9.4-Python-2.7.13/data/insecta_odb9 -m genome -f



