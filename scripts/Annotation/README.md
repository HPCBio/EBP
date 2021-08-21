# 1. Background

Gene prediction or gene finding refers to the process of identifying the regions of genomic DNA that encode genes. This includes protein-coding genes as well as RNA genes, but may also include prediction of other functional elements such as regulatory regions. Gene finding is one of the first and most important steps in understanding the genome of a species once it has been sequenced. 

<p>
<img align="left" src="./docs/gene_prediction_explained.jpg" />

</br></br></br>
</p>

Source: https://en.wikipedia.org/wiki/Gene_prediction


We perform gene prediction with the  GeneMark-EP + AUGUSTUS + BRAKER2 pipeline. This pipeline combines extrinsic and ab initio approaches by mapping protein and EST data to the genome to validate ab initio predictions. Augustus results can also incorporate hints in the form of EST alignments or protein profiles to increase the accuracy of the gene prediction. 

See: http://opal.biology.gatech.edu/GeneMark/

Then we perform protein function prediction using InterProScan and blastp with the UniprotDB/swissprot.

See: 

https://interpro-documentation.readthedocs.io/en/latest/interproscan.html

https://www.uniprot.org/


# 2. Dependencies

- This step should be run after purge_dups, scaffolding, gap filling and masking

- Orthologs database from https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz

- ProtHint/2.5.0 

- key files provided with this repo: .gm_key .gm_key_64

- BRAKER 2.1.6

- cdbfasta 20181005

- Cufflinks 2.2.1

- InterProScan 5.47-82.0

- BLAST+ 2.10.1

- Uniprot and swissprot databases in fasta and blast-ready formats. https://www.uniprot.org/downloads#uniprotkblink

- Perl

- scritps provided with this repo: maker_functional_fasta, maker_functional_gff, ipr_update_gff, BUSCO_protein_slurm.sh 

- BUSCO

- assemblathon

- slurm workload manager https://slurm.schedmd.com/documentation.html


# 3. Gene Prediction

The script prothint_genemark_braker_slurm.sh runs the gene prediction pipeline.

## 3.1 Setup

- Copy the key files to your home folder. GeneMark looks for them in the $PATH global variable

- Download a curated database of close orthologs in fasta format. 
We chose https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz
Must edit the database and replace ambiguous amino acids such as * with Ns.
The gene predictions will be as good as this databasel so it is essential 
that you select a set of curated orthologs that is closest to your assembly in the taxonomy tree.

- Sequence headers with long names should be avoided. 
Make sure that your database of orthologs and your genome have short headers.

## 3.2 Run the gene prediction script

Edit the script prothint_genemark_braker_slurm.sh with your information

To run the script

<pre>
sbatch prothint_genemark_braker_slurm.sh
</pre>

The execution could take from 4 hours to more than a day.

## 3.3 Outputs

<pre>
braker/
cufflinks/
EB15_shortHeaders.fa.gz
genemark/
prothint/
</pre>


- Gene models are found in prothint and genemark folders

- proteome in fasta and gff formats are in cufflinks folder


## 3.4 run QC 

The script BUSCO_protein.sh generates proteome completeness report calculated with BUSCO.

# 4. Protein Function Prediction

We used InterProScan and blastp with the UniprotDB/swissprot databases

## 4.1. Setup

The script is called protein_annotation_slurm.sh and it requires other programs that are provided with this repo.

- Copy the scripts to your path

- Index the databases for blast.

- Edit the script with your information.

- The input file is called proteins.fa that the previous step produced. It is located inside the cufflinks folder.

## 4.2 Run the script


<pre>
sbatch protein_annotation_slurm.sh
</pre>

The execution could take a few hours.

## 4.3 Output

Example:


<pre>
braker/
BUSCO_proteins_vs_insecta/
cufflinks/
EB15_output.blastp
EB15_output.iprscan
EB15_proteins.putative_function.domain_added.gff
EB15_proteins.putative_function.fasta.gz
EB15_proteins.putative_function.gff
EB15_proteins.putative_function_rename.gff
EB15_shortHeaders.fa.gz
genemark/
prothint/
</pre>


The BUSCO results are inside the  BUSCO_proteins_vs_insecta/ folder

The annotated proteome in fasta and gff format are the files: 

- EB15_proteins.putative_function.fasta.gz

- EB15_proteins.putative_function_rename.gff






