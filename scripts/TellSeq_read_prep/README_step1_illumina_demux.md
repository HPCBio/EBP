# Background on TellSeq libraries 

TellSeq libraries: https://www.universalsequencing.com/technology


TellSeq libraries produce short reads.  These short reads were used in this project for two different purposes: 

1. To estimate genome size using kmer distribution histograms of the short reads and then running GenomeScope 2.0 http://qb.cshl.edu/genomescope/genomescope2.0/

2. To scaffold the assembly with short reads.  The scaffolding toolkit we used was ARCS-LINKS. See page: https://github.com/bcgsc/arcs


# Dependencies

This program expects the following tools/languages to be installed as modules and be available in your path:

- bcl2fastq  version 2.20

# Inputs

The inputs are scattered in separate files and folders inside the rundir that the Illumina instrument created to perform the sequencing task.

For more details see bcl2fastq v2.20 Software Guide.



## Samplesheet

The samplesheet is a text file describing the samples to be demultiplexed.

For more details see  https://support.illumina.com/help/BaseSpace_Sequence_Hub/Source/Informatics/BS/SampleSheets_swBS.htm

<pre>

$ tail /rundir/Data/Intensities/BaseCalls/SampleSheet_TellSeq.csv

[Settings],,,,,,,,,,
ReverseComplement,0,,,,,,,,,
Adapter,,,,,,,,,,
AdapterRead2,,,,,,,,,,
,,,,,,,,,,
[Data],,,,,,,,,,
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,User_Name,Report_No,Application,Primers,Primer_Mismatch
3,Sample_EB5_Tellseq,EB5_Tellseq,,AAGGTTCA,Project_RobTellSeq,robinson,1,,,
3,Sample_EB14_Tellseq,EB14_Tellseq,,ACTTAGCA,Project_RobTellSeq,,,,,
3,Sample_EB20_Tellseq,EB20_Tellseq,,TGTTCTAG,Project_RobTellSeq,,,,,

</pre>

## Demultiplexing

The demultiplexing command has to be launched from /rundir/Data/Intensities/BaseCalls/

The param --tiles s_3 is specified in the cmd to indicate that only lane 3 will be demultiplexed with this command

<pre>

$ cd /rundir/Data/Intensities/BaseCalls/

$ nohup bcl2fastq -R /rundir/  \
--processing-threads 24  \
--output-dir /rundir/TellSeq_notrim \
--sample-sheet /rundir/Data/Intensities/BaseCalls/SampleSheet_TellSeq.csv \
--create-fastq-for-index-reads \
--barcode-mismatches 0 \
--ignore-missing-bcls \
--use-bases-mask Y150n,Y18,I8n*,Y150n \
--minimum-trimmed-read-length 1 \
--mask-short-adapter-reads 1 \
--no-bgzf-compression \
--tiles s_3  > nohup_TellSeq.log &

</pre>

# Outputs

The output folder will have as many projects and samples as described in the sample sheet.

For example; the samplesheet shown earlier would produce an output folder with this structure

<pre>

$ tree
.
|-- Sample_EB14_Tellseq
|   |-- EB14_Tellseq_S3_L003_I1_001.fastq.gz
|   |-- EB14_Tellseq_S3_L003_R1_001.fastq.gz
|   |-- EB14_Tellseq_S3_L003_R2_001.fastq.gz
|   |-- EB14_Tellseq_S3_L003_R3_001.fastq.gz
|-- Sample_EB20_Tellseq
|   |-- EB20_Tellseq_S4_L003_I1_001.fastq.gz
|   |-- EB20_Tellseq_S4_L003_R1_001.fastq.gz
|   |-- EB20_Tellseq_S4_L003_R2_001.fastq.gz
|   |-- EB20_Tellseq_S4_L003_R3_001.fastq.gz
|-- Sample_EB5_Tellseq
|   |-- EB5_Tellseq_S5_L003_I1_001.fastq.gz
|   |-- EB5_Tellseq_S5_L003_R1_001.fastq.gz
|   |-- EB5_Tellseq_S5_L003_R2_001.fastq.gz
|   `-- EB5_Tellseq_S5_L003_R3_001.fastq.gz

</pre>

Each sample will have a separate folder and should contain these FOUR files

R1 has forward reads, readlen=150

R3 has reverse reads, readlen=150

R2 has 10x-like barcodes, readlen=18

I1 has illumina indices, readlen=8. This file is not used during post-processing.
