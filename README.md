# iRibo

A comprehensive tool for ... (brief description of what iRibo does).

------------------------------------------------------------------------------

## Table of Contents:
1. Installation and Prerequisites
2. Data Collection and Preprocessing
3. GetCandidateORFs
4. GenerateTranslationProfile
5. GenerateTranslatome

------------------------------------------------------------------------------

## Installation and Prerequisites:

- iRibo Download: You can obtain iRibo from here: https://github.com/Alistair-Turcan/iRibo. 
  Follow these steps after downloading:
  1. Download the zip file.
  2. cd into the directory.
  3. Run make.
  4. Note: It is compiled using the g++ compiler and c++17.

- System Requirements:
  - Memory: At least 64GB
  - Storage: At least 100GB

- Dependencies: Ensure you have R version 4.2.2 installed. 
  You can download and access the documentation here: https://www.r-project.org/.

------------------------------------------------------------------------------

## Data Collection and Preprocessing:

iRibo requires the following inputs:
1. Genome: In FASTA format.
2. Annotations: Accepts GFF3 or GTF formats.
3. Transcriptome (optional): Contains both annotated and unannotated transcripts.
4. Aligned ribo-seq reads: Either in SAM or BAM format.

Recommendation: 
- Ribo-seq samples in fastq format should have low-quality reads and adaptors trimmed for optimal alignment.
- Ensure the genome, annotations, and transcriptome have matching chromosome identifiers, e.g., >chr1 in the genome matches chr1 in the annotations.

------------------------------------------------------------------------------

## GetCandidateORFs:

Start with creating candidate ORFs for assessing translation. The result is a file named candidate_orfs containing potential ORFs for iRibo's evaluation.

Options:
- --Transcriptome=path/to/transcriptome.gtf: Use transcriptome instead of the genome for creating ORFs.
- --Output=path/to/output_folder: Designate the output directory.
- --Threads=1: Define the number of threads.

------------------------------------------------------------------------------

## GenerateTranslationProfile:

Generate a genome-wide profile of translation with aligned ribo-seq reads. This stage produces several files:
- translation_calls: Statistics on reads per ORF.
- null_distribution: Null distribution stats of reads per ORF.
... (list the rest of the files in the same manner)

SAMs Construction:
Construct sams.txt that holds a list of paths to all SAM/BAM files, separated by line:
sam_dir/SRR1042853_aligned.out.bam
sam_dir/SRR1042855_aligned.out.bam
...

Options:
- --Output=path/to/output_folder: Set output directory.
- --Threads=1: Number of threads.
... (list the rest of the options in the same manner)

------------------------------------------------------------------------------

## GenerateTranslatome:

This final step produces the translatome, resulting in two files:
- translated_orfs.csv: Data on detected ORFs as translated.
- nORF_discovery.png: Graphical representation of p-values for actual and scrambled ORFs, alongside the FDR cutoff.

Options:
- --Output=path/to/output_folder: Define output directory.
- --Threads=1: Determine the number of threads.
... (list the rest of the options in the same manner)

------------------------------------------------------------------------------

You can add a Contribution, License, or any other general section as required.

------------------------------------------------------------------------------
