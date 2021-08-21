# Background on TellSeq libraries 

TellSeq libraries: https://www.universalsequencing.com/technology


TellSeq libraries  were used in this project for two different purposes: 

1. To estimate genome size using kmer distribution histograms of the short reads and then running GenomeScope 2.0 http://qb.cshl.edu/genomescope/genomescope2.0/

2. To scaffold the assembly with short reads.  The scaffolding toolkit we used was ARCS-LINKS. See page: https://github.com/bcgsc/arcs


## Read preparation steps

- [Step1: demultiplex raw bcl files to fastq files](README_step1_illumina_demux.md)

- [Step2: Install TellSeq container on cluster](README_step2_TellSeq_Installation_on_biotech.md)

- [Step3: Run TellRead to QC-trim reads](README_step3_run_Tellread.md)

## Reformatting steps

- [Step4: Reformat TellSeq reads for supernova toolkit](README_step4_reformat4supernova.md)

- [Step5: Reformat TellSeq reads for ARCS](README_step5_reformat4arcs.md)

## Genomescope 

- [Step6: Using TellSeq reads to estimate genome size with GenomeScope ](README_step6_run_GenomeScope.md)




