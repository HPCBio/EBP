#!/bin/bash

# ----------------SLURM Parameters----------------
#SBATCH -J 20210531_EB20_FreeBayes_align_to_consensus #jobname
#SBATCH --mem=100g
#SBATCH -N 1 #nodes
#SBATCH -n 12 #cpus-per-task
#SBATCH -p normal #queue
#SBATCH -A ebp
#SBATCH --mail-user=kosterbu@illinois.edu
#SBATCH --mail-type=ALL
#SBATCH -D /home/groups/earthbiogenome/results/20210524_EB20_FreeBayes/

# ----------------Commands------------------------

# Link files to working directory:

#ln -s /home/groups/earthbiogenome/results/20210428_EB20_RAILS_SEALER/EB20_Salsa_Sealer_scaffold.fa .

#ln -s /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz .
#ln -s /home/groups/earthbiogenome/data/TellSeq_not_interleaved_reads/EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz .

# ----------------Load Modules--------------------

module load  BWA/0.7.17-IGB-gcc-8.2.0 

bwa index  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa

bwa mem -t 12 EB20_Salsa_Sealer_scaffold_rename_cleaned.fa EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz >  EB20_purged_asm_vs_TellSeq.sam

# ----------------Load Modules--------------------

module purge

module load SAMtools/1.12-IGB-gcc-8.2.0

# ----------------Commands------------------------


samtools view -@ 12 -o  EB20_purged_asm_vs_TellSeq.bam  EB20_purged_asm_vs_TellSeq.sam

samtools sort -@ 12 -o  EB20_purged_asm_vs_TellSeq_sorted.bam  EB20_purged_asm_vs_TellSeq.bam

samtools index  EB20_purged_asm_vs_TellSeq_sorted.bam

samtools faidx  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa

# ----------------Load Modules--------------------

module purge
module load Python/3.7.2-IGB-gcc-8.2.0
module load FreeBayes/1.3.4-IGB-gcc-8.2.0

# ----------------Commands------------------------

# FreeBayes examples:

freebayes -f  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa  EB20_purged_asm_vs_TellSeq_sorted.bam >  EB20_purged_asm_vs_TellSeq_sorted.vcf

# ----------------Load Modules--------------------
module purge
module load VCFtools/0.1.16-IGB-gcc-8.2.0-Perl-5.28.1 
module load BCFtools/1.12-IGB-gcc-8.2.0
module load Merfin/20210507-IGB-gcc-8.2.0
module load meryl/1.3
module load Jellyfish/2.3.0-IGB-gcc-8.2.0

#vcf=P350  # keep short without dots (.)
#genome=P350_chromosomal.fasta
#reads=P350R*.fastq
#readR1=P350R1.fastq
#readR2=P350R2.fastq
#kmer=21
#ploidy=2

## Normalize VCF file output from freebayes and convert to bcf
bcftools --version
bcftools view -Ou -e'type="ref"'  EB20_purged_asm_vs_TellSeq_sorted.vcf | bcftools norm -Ob -f  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa -o  EB20_purged_asm_vs_TellSeq_norm.bcf --threads 12

## Filter normalized bcf file and output as vcf

bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=12  EB20_purged_asm_vs_TellSeq_norm.bcf >  EB20_purged_asm_vs_TellSeq_filtered.vcf.gz
bcftools sort -o  EB20_purged_asm_vs_TellSeq_filtered_sorted.vcf.gz -m 10G -Oz  EB20_purged_asm_vs_TellSeq_filtered.vcf.gz
bcftools index  EB20_purged_asm_vs_TellSeq_filtered_sorted.vcf.gz


## create Meryl databases for input to merfin

meryl count threads=12 k=19 EB20_Salsa_Sealer_scaffold_rename_cleaned.fa  output  EB20_purged.meryl
#meryl count k=19 EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz output EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.meryl
#meryl count k=19 EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz output EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.meryl
#meryl union-sum EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.meryl EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.meryl output EB20_nobc_all.meryl

## Run JellyFish to get histogram peak as input to merfin

#jellyfish --version
# run jellyfish to get kmer data
#jellyfish count -C -m 21 -s 3000000000 -t 30 <(zcat EB20_nobc_R1_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz) <(zcat EB20_nobc_R2_T500.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz) -o EB20_reads.jf
# create histogram from that data
#jellyfish histo -t 36 EB20_reads.jf > EB20_reads.hist
# run genomescope with hist

#genomescope.R -i ${vcf}reads.hist -k 21 -o P348_genomescope_out -p 2 --fitted_hist

# identify peak
head -n 50 EB20_reads.hist
peak=`more +5 EB20_reads.hist | sort -k 2n | tail -n 1 | awk '{print $1}'`


## Merfin filter with meryl counts and jellyfish peak

merfin -threads 12 -vmer -sequence  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa -seqmers  EB20_purged.meryl -readmers EB20_nobc_all.meryl -peak ${peak} \
-vcf  EB20_purged_asm_vs_TellSeq_filtered_sorted.vcf.gz  \
-output  EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz >  EB20_purged_asm_vs_TellSeq_filtered_sorted_merfin.out 2>&1


## bcftools build consensus using filtered vcf file

bcftools view -Oz --threads=12 EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz.polish.vcf >  EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz.polish.vcf.gz #bgzip merfin output

bcftools index  EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz.polish.vcf.gz

bcftools consensus  EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz.polish.vcf.gz -f  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa -H 1 >  EB20_purged_FreeBayes_merfin_polished.fasta 
# -H 1 applies only first allele from GT at each positionscrips

## Report number of changes:

bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Ov EB20_purged_asm_vs_TellSeq_filtered_sorted_out.dump.gz.polish.vcf.gz | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > EB20.numvar
echo "Num. bases affected: `cat EB20.numvar`"

perl /home/a-m/kosterbu/scripts_misc/get_fasta_stats_102010.pl -g -T   EB20_purged_FreeBayes_merfin_polished.fasta
perl /home/a-m/kosterbu/scripts_misc/N50.pl   EB20_purged_FreeBayes_merfin_polished.fasta

#---------------------VGP method for consensus---------------

bcftools consensus -Hla -f  EB20_Salsa_Sealer_scaffold_rename_cleaned.fa  EB20_purged_asm_vs_TellSeq_filtered_sorted.vcf.gz >  EB20_purged_FreeBayes_VGP_polished.fasta

perl /home/a-m/kosterbu/scripts_misc/get_fasta_stats_102010.pl -g -T  EB20_purged_FreeBayes_VGP_polished.fasta
perl /home/a-m/kosterbu/scripts_misc/N50.pl  EB20_purged_FreeBayes_VGP_polished.fasta
