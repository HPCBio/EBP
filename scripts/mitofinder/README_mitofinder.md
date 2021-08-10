# 1. Background

We used MitoFinder to identify and annotate mtDNA.  

git repo: https://github.com/RemiAllio/MitoFinder

MitoFinder: Efficient automated large-scale extraction of mitogenomic data in target enrichment phylogenomics
Mol Ecol Resour. 2020 Jul; 20(4): 892–905.
Published online 2020 Apr 25. doi: 10.1111/1755-0998.13160

# 2. Dependencies

- singularity

- MitoFinder v 1.4

- reference mitochondrial genome in GenBank format  https://www.ncbi.nlm.nih.gov/genome/browse#!/organelles/

- slurm workload manager https://slurm.schedmd.com/documentation.html

# 3. Setup

MitoFinder needs not only the target genome file, but also the reference mitochondrial genome in GenBank format.
This step cannot be automated. It is basically a search for the mito genome that is closest to the target genome.
We only provide the NCBI download page for all mitochondrial genomes they have.

# 4. Run the script

MitFinder is a container that can be run with singularity.

Edit the script with your information before launching it.

# 5. Check the results

## 5.1 Check the log

MitoFinder writes a log file with a summary of the analysis.

Example:

<pre>
...
MitoFinder found a single mitochondrial contig
Checking resulting contig for circularization...

Evidences of circularization were found!
Sequence is going to be trimmed according to circularization position. 



Creating summary statistics for the mtDNA contig

Annotating

tRNA annotation with MitFi run well.

Annotation completed


Creating GFF and fasta files.

Note: 
15 genes were found in mtDNA_contig

</pre>

## 5.2 Inspect the output files. 

Example:

<pre>
EB1_MitoFinder_whole_genome_final_genes_AA.fasta.gz
EB1_MitoFinder_whole_genome_final_genes_NT.fasta.gz
EB1_MitoFinder_whole_genome.infos
EB1_MitoFinder_whole_genome_mtDNA_contig.fasta.gz
EB1_MitoFinder_whole_genome_mtDNA_contig.gb
EB1_MitoFinder_whole_genome_mtDNA_contig_genes_AA.fasta.gz
EB1_MitoFinder_whole_genome_mtDNA_contig_genes_NT.fasta.gz
EB1_MitoFinder_whole_genome_mtDNA_contig.gff
EB1_MitoFinder_whole_genome_mtDNA_contig.tbl

</pre>


# 6. DNA barcoding

<pre>
The BOLD Identification System (IDS) for COI accepts sequences from 
the 5' region of the mitochondrial Cytochrome c oxidase subunit I gene 
and returns a species-level identification when one is possible. 
Further validation with independent genetic markers will be desirable in some forensic applications. 
</pre>

https://www.boldsystems.org/index.php/IDS_OpenIdEngine

Include here how the COI fasta sequence is recovered from the MitoFinder results

# 7. Visualization

Load the file to  http://wolfe.ucd.ie/GenomeVx/

Include here which file is loaded to this web server





