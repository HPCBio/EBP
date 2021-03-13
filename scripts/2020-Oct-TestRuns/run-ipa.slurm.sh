#!/bin/bash
#SBATCH -n 49
#SBATCH --mem=500g
#SBATCH -p normal 
#SBATCH --mail-type=END,FAIL
#SBATCH -J pbipa

### Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

### Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

### Load Modules
module load pbipa/1.3.0

### Run app on file

ipa local --nthreads 12 \
     --njobs 4 \
     -i $FASTQ \
     --resume
