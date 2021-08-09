# 1. Background

purge_dups purges/deletes haplotigs and overlaps in an assembly based on read depth.

See https://github.com/dfguan/purge_dups

# 2. Dependencies

- This step should be run after a raw assembly is generated with HiFi reads

- purge_dups 1.2.5

- slurm workload manager https://slurm.schedmd.com/documentation.html

# 3. Setup

## 3.1 create a text file in the working directory called "pbfofn" listing the complete path(s) of the PacBio fastq file, one per line

<pre>
cd workingdir

ln -s /path/to/hifi_1_consensus.fastq
ln -s /path/to/hifi_2_consensus.fastq

ls -1 *fastq > pbfofn

</pre>

## 3.2  run pd_config.py to generate a config file for the purge_dups run, make sure you have correct path to reference fasta file

<pre>

usage: pd_config.py [-h] [-s SRF] [-l LOCD] [-n FN] [--version] ref pbfofn

generate a configuration file in json format

positional arguments:
  ref                   reference file in fasta/fasta.gz format
  pbfofn                list of pacbio file in fastq/fasta/fastq.gz/fasta.gz format (one absolute file path per line)

optional arguments:
  -h, --help            show this help message and exit
  -s SRF, --srfofn SRF  list of short reads files in fastq/fastq.gz format (one record per line, the
                        record is a tab splitted line of abosulte file path
                        plus trimmed bases, refer to
                        https://github.com/dfguan/KMC) [NONE]
  -l LOCD, --localdir LOCD
                        local directory to keep the reference and lists of the
                        pacbio, short reads files [.]
  -n FN, --name FN      output config file name [config.json]
  --version             show program's version number and exit
  
</pre>

This step should create a config.json file as specified with the --name parameter


## 3.3 Edit config.json

Open the config.json file and edit lines to skip BUSCO and KCM steps.

In other words change busco and kcm blocks FROM "skip": 0  TO "skip": 1

# 4. Run the script

The script that performs this step is called purge_dups.slurm.sh

Edit it according to the information of your genome.

To run the script

<pre>
sbatch purge_dups.slurm.sh
</pre>

# 5. Outputs

After the pipeline is finished, there will be four new directories in the working directory (set in the configuration file).

- coverage: coverage cutoffs, coverage histogram and base-level coverage files

- split_aln: segmented assembly file and a self-alignment paf file.
    
- purge_dups: duplicate sequence list.

- seqs: purged primary contigs ending with .purge.fa and haplotigs ending with .red.fa, also K-mer comparison plot and busco results are also in this directory.

## 6. QC the assembly

We provide a script BUSCO_metrics_slurm.sh that runs assemblathon.pl to calculate assembly contiguity.

It also runs BUSCO (database: insecta_odb10)  to calculate genome completeness.




