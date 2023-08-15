# iRibo

## Installation and Prerequisites

iRibo can be downloaded from https://github.com/Alistair-Turcan/iRibo. Simply download the zip file, cd into the directory, and run make. It is compiled using the g++ compiler and c++17.

A computer with at least 64GB of memory and 100GB of storage is typically required.

R version 4.2.2 is required, Download and documentation is available at https://www.r-project.org/.

Data Collection and Preprocessing 

iRibo has 4 inputs:

Genome, in FASTA format.
Annotations, in either GFF3 or GTF format.
Transcriptome (optional). Contains all transcripts, both annotated and unannotated.
Aligned ribo-seq reads in either SAM or BAM format.
Ribo-seq samples in fastq format should have low quality reads and adaptors trimmed for optimal alignment.

Make sure the genome, annotations, and transcriptome all have identical chromosome identifiers. Ie, >chr1 in the genome, and chr1 in the annotations.

## GetCandidateORFs

The first step of iRibo is to create candidate ORFs to assess translation in. This step will create a file, candidate_orfs, containing all the potential ORFs that iRibo will look for translation patterns in.


Additional options include:
--Transcriptome=path/to/transcriptome.gtf. Creates candidate ORFs from the transcriptome instead of the genome.
--Output=path/to/output_folder. Specifies output directory.
--Threads=1. Specifies number of threads to use.



## GenerateTranslationProfile

The next step is to create a genome-wide profile of translation using the aligned ribo-seq reads. This step will create a number of files:
translation_calls contains statistics about the reads in each orf
null_distribution contains statistics about a null distribution of the reads in each orfs
all_passed_reads_f contains all reads that passed quality control on the forward strand, whether or not they are inside a candidate orf.
all_passed_reads_r contains the same as all_passed_reads_f, but for the reverse strand.
candidate_orfs.gff3 contains annotations of all ORFs.
all_passed_reads_f.wig contains tracks of all the forward strand reads.
all_passed_reads_r.wig contains tracks of all the reverse strand reads.



sams.txt is a user-constructed file that contains a list of filepaths to all of the SAM/BAM files, separated by line. It should look like:

sam_dir/SRR1042853_aligned.out.bam
sam_dir/SRR1042855_aligned.out.bam
sam_dir/SRR1042857_aligned.out.bam
…

Additional options for this step include:
--Output=path/to/output_folder. Specifies output directory.
--Threads=1. Specifies number of threads to use.
--Min_Length=25. Minimum read length to perform quality control in. All reads below this length will be discarded.
--Max_Length=35. Maximum read length to perform quality control in. All reads above this length will be discarded.
--P_Site_Distance=20. Max distance to check for a p-site. Most are below 20.
--QC_Count=10000. How many reads must be in the first frame of canonical genes in a read length in order to pass quality control.
--QC_Periodicity=2.0. The scale of periodicity in canonical genes required in a read length in order to pass quality control.
--QC_Positions=false. Whether to use positions or read counts in quality control. Positions is better for sparser data, read counts is better for more dense data.

## GenerateTranslatome

The last step is to generate the translatome. This produces two files. translated_orfs.csv, containing information about all the ORFs iRibo detects as translated, and nORF_discovery.png, showing a graph of the p-values of the actual and scrambled ORFs, as well as the FDR cutoff.



Additional options for this step include:
--Output=path/to/output_folder. Specifies output directory.
--Threads=1. Specifies number of threads to use.
--ExcludeChr=chr1,chr8. Specifies which chromosome/contigs to exclude by a comma separated list of id’s, if any.
--ExcludeOverlapGene=True. Specifies to exclude nORFs that overlap canonical genes on the same strand.
--FDR=0.05. Specifies desired false discovery rate.
