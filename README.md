------------------------------------------------------------------------------

# iRibo

A comprehensive tool for integrating ribosome profiling data to detect genome wide translation.

------------------------------------------------------------------------------

## Table of Contents:
1. Installation and Prerequisites
2. Data Collection and Preprocessing
3. GetCandidateORFs
4. GenerateTranslationProfile
5. GenerateTranslatome

------------------------------------------------------------------------------

## Installation and Prerequisites:

- iRibo Download: 
  Obtain iRibo from https://github.com/Alistair-Turcan/iRibo. 
  Instructions post download:
  1. Download the zip file.
  2. cd into the directory.
  3. Run make.
  Note: Compiled using g++ compiler and c++17.

- System Requirements:
  - Memory: At least 64GB
  - Storage: At least 100GB

- Dependencies: 
  R version 4.2.2 is necessary. Download and documentation are available at https://www.r-project.org/.

------------------------------------------------------------------------------

## Data Collection and Preprocessing:

iRibo takes in:
1. Genome in FASTA format.
2. Annotations in GFF3 or GTF format.
3. Transcriptome (optional): Accepts both annotated and unannotated transcripts in GFF3 or GTF format.
4. Aligned ribo-seq reads in SAM or BAM format.

Recommendation: 
- Trim low-quality reads and adapters from ribo-seq samples in FASTQ format for best alignment.
- Choose a transcriptome that is as comprehensive as possible, for both read mapping and iRibo.

Note:
- Ensure genome, annotations, and transcriptome have consistent chromosome identifiers, e.g., >chr1 in genome should match chr1 in annotations.

------------------------------------------------------------------------------

## GetCandidateORFs:

Begin by generating candidate ORFs for translation assessment. This results in three files:
- candidate_orfs: A list of all candidate ORFs that will be assessed for translation. This is filtered due to overlapping translation signals being hard to separate.
- all_orfs: A list of every possible ORF that exists in the genome or transcriptome.
- candidate_orfs.gff3: Annotations of all candidate ORFs, ready to be put into a genome browser like IGV.

./iRibo --RunMode=GetCandidateORFs --Genome=path/to/genome.fa --Annotations=path/to/annotations.gtf

Options:
- --Transcriptome=path/to/transcriptome.gtf: Use transcriptome over genome.
- --Output=path/to/output_folder: Specify output directory.
- --Threads=1: Set number of threads.

The columns of candidate_orfs and all_orfs can be interpreted as such:
- CandidateORF_ID: index of the ORF
- Transcript_ID: ID of the transcript the ORF is on
- Gene_ID: The name of the gene, if canonical. X otherwise
- contig: Index form of contig/chromosome the ORF is on.
- strand: Strand the ORF is on.
- ORF_coord1: First coordinate of the ORF
- ORF_coord2: Second coordinate of the ORF
- genomic_coordinates: Coordinates of the ORF. The start and stop position within each exon of the ORF, separated by -.
- ORF_length: Exonic length of the ORF.
- antisense_gene: Protein-coding gene on the antisense strand, if any. X otherwise.
- gene_intersect: Protein-coding gene overlap, if any. X otherwise.
- contig_str: String form of contig/chromosome the ORF is on.


------------------------------------------------------------------------------

## GenerateTranslationProfile:

Generate a genome-wide translation profile using aligned ribo-seq reads. This phase outputs:
- translation_calls: Read statistics per ORF.
- null_distribution: Null distribution read statistics per ORF.
- all_passed_reads_f: Quality-passed forward strand reads.
- all_passed_reads_r: Quality-passed reverse strand reads (similar to forward).
- all_passed_reads_f.wig: Tracks of forward strand reads.
- all_passed_reads_r.wig: Tracks of reverse strand reads.

To run:
./iRibo --RunMode=GenerateTranslationProfile --Genome=path/to/genome.fa --Riboseq=path/to/sams.txt --CandidateORFs=path/to/candidate_orfs

SAMs File:
sams.txt should list paths to all SAM/BAM files, separated by lines:

e.g., 

sam_dir/SRR1042853_aligned.out.bam

sam_dir/SRR1042855_aligned.out.bam

Options:
- --Output=path/to/output_folder: Define output directory.
- --Threads=1: Set thread count.
- --Min_Length=25: Minimum read length for quality control. Most riboseq reads are between length 25-35nt, but can vary by experiment.
- --Max_Length=35: Maximum read length for quality control. Most riboseq reads are between length 25-35nt, but can vary by experiment.
- --P_Site_Distance=20: Max distance to check for a P-site. Most p-sites are below 20. If you have a reason to think it could be longer, can be increased.
- --QC_Count=10000: Number of reads in the first frame of protein-coding genes for quality control.
- --QC_Periodicity=2.0: Periodicity scale in canonical genes for quality control.
- --QC_Positions=false: Use positions or read counts in quality control.

------------------------------------------------------------------------------

## GenerateTranslatome:

This final step creates the translatome, yielding:
- translated_orfs.csv: Data on translated ORFs. Adds 3 new columns for in-frame read count, expression level (length/in-frame read count), and p-value of discovery.
- nORF_discovery.png: Graph for p-values of real vs. scrambled nORFs and FDR cutoff.
- cORF_discovery.png: Graph for p-values of real vs. scrambled cORFs and FDR cutoff.
- translated_orfs.gff3: Annotations of all translated ORFs, ready to be put into a genome browser like IGV.

To run:
Rscript GenerateTranslatome.R --TranslationCalls=path/to/translation_calls --NullDistribution=path/to/null_distribution --CandidateORFs=path/to/candidate_orfs

Options:
- --Output=path/to/output_folder: Designate output directory.
- --Threads=1: Specify thread count.
- --ExcludeChr=chr1,chr8: Chromosomes/contigs to exclude. Default is none.
- --ExcludeOverlapGene=True: Exclude nORFs overlapping canonical genes on the same strand.
- --FDR=0.05: Define desired false discovery rate.
- --Scrambles=100: Set number of scrambles to calculate FDR.

------------------------------------------------------------------------------
