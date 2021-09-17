# Earth Biogenome Project (EBP) Insect Genome Assembly


The Illinois Innovation Network and the Discovery Partners Institute funded a pilot project to assemble the genomes of agriculturally relevant insects in Illinois for which little or no genomic data are available. Hiqh-quality DNA was isolated using protocols optimized for small, difficult samples. The pilot project Implemented the novel use of Tell-Seq linked-read libraries for the dual purpose of genome size estimation and linked-read scaffolding. Genomes were assembled using PacBio HiFi reads, Tell-Seq reads, and Dovetail Omni-C reads for chromosome-range scaffolding. 

Eight high-quality genomes were assembled from non-model organisms with contig N50 >1Mb and scaffold N50 >5Mb, including the second-only soon-to-be public genome for the order Neuroptera. Species were confidently identified using the mitochondrial genomes assembled from the HiFi reads. Potential endosymbionts and pathogens were identified as well as novel prey information from predator species. Genomes were annotated using the BRAKER2 pipeline, generating a rich set of novel data to mine

Our sponsors: 

- Illinois Innovation Network https://iin.uillinois.edu/

- Discovery Partners Institute https://dpi.uillinois.edu/

- Carl R. Woese Institute for Genomic Biology (IGB) https://www.igb.illinois.edu

- Roy J. Carver Biotechnology Center https://biotech.illinois.edu/


To learn more about the Earth Biogenome Project, please use this link: https://www.earthbiogenome.org/


# The Workflow

<p>
<img align="left" src="./docs/EBP_Workflow_1.png" />

</br></br></br>
</p>

<p>
</br></br></br>
</p>



# Denovo genome assembly using HiFi reads

These are the steps:


1. [Generate raw assembly with hifiasm](/scripts/Raw_assembly/README.md)

2. [Purge duplicate contigs](/scripts/purge_dups/README.md)

3. [Scaffold using TellSeq reads](/scripts/scaffolding/README.md)

4. [Scaffold using Omni-C reads](/scripts/scaffolding/README.md)

5. [Fill gaps](/scripts/gap_filling_and_masking//README.md)

6. [Mask repeats and low complexity regions](/scripts/gap_filling_and_masking/README.md)

7. [Predict gene and protein function](/scripts/Annotation/README.md)

8. [Identify and annotate mitochondrial DNA](/scripts/mitofinder/README.md)

9. [Assess genome completeness w Merqury](/scripts/Merqury_completeness/README.md)

10. [Identify contaminants and artifacts in genome](/scripts/blobtools_contaminants_detection/README.md)



# Denovo genome assembly using CLR reads

1. [Generate raw assembly with Redbean](/scripts/Raw_assembly//README.md)

2. [Base-correct assembly with Arrow](/scripts/Arrow_polish/README.md)

3. [Purge duplicate contigs](/scripts/purge_dups/README.md)

4. [Pilon polishing](/scripts/pilon_polishing/README.md)

5. [Scaffold using TellSeq reads](/scripts/scaffolding/README.md)

6. [Scaffold using Omni-C reads](/scripts/scaffolding/README.md)

7. [Fill gaps](/scripts/gap_filling_and_masking//README.md)

8. [Mask repeats and low complexity regions](/scripts/gap_filling_and_masking/README.md)

9. [Predict gene and protein function](/scripts/Annotation/README.md)

10. [FreeBayes polishing](/scripts/FreeBayes_polishing/README.md)

11. [Identify and annotate mitochondrial DNA](/scripts/mitofinder/README.md)

12. [Assess genome completeness w Merqury](/scripts/Merqury_completeness/README.md)

13. [Identify contaminants and artifacts in genome](/scripts/blobtools_contaminants_detection/README.md)

