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
  R version 4.2.2 or above is necessary. Download and documentation are available at https://www.r-project.org/.

  Samtools. Download and documentation are available at https://www.htslib.org/. 

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
- riboseq_reads_plus.wig: Tracks of forward strand reads.
- riboseq_reads_minus.wig: Tracks of reverse strand reads.

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

The columns of translation_calls are:
- index: The index (0-based) of the orf in the candidate_orfs file used.
- frame0: How often the first frame is the biggest in the orf, when any frame is the biggest. This is the number of successes in the binomial test for translation.
- frame_sum: The total sum of how often any frame is the biggest in the orf, if any. This is the number of trials in the binomial test for translation.
- reads0, reads1, reads2: The total count of reads mapping to the first, second, and third frames of the orf, respectively.

The columns of null_distribution are:
- scrambledX: Where X is a number 0-99, and contains how often the first frame is the biggest in that scramble of the orf, when any frame is the biggest.
- scrambled_sumX: The total sum of how often any frame is the biggest in that scrambled orf, if any.

The columns of all_passed_reads_f and all_passed_reads_r are:
- chr: The index of the chromosome for the read, in the order that chromosome appears in the genome.
- strand: The strand the read aligns to, 0 for forward, 1 for reverse.
- pos: The genomic position of the start of the read.
- count: How many reads map to this chr, strand, pos.


------------------------------------------------------------------------------

## GenerateTranslatome:

This final step creates the translatome, yielding:
- translated_orfs.csv: Data on translated ORFs. The same output format as candidate_orfs, but with 3 new columns for in-frame read count, expression level (length/in-frame read count), and p-value of discovery.
- nORF_discovery.png: Graph for p-values of real vs. scrambled nORFs and FDR cutoff.
- cORF_discovery.png: Graph for p-values of real vs. scrambled cORFs and FDR cutoff.
- translated_orfs.gff3: Annotations of all translated ORFs, ready to be put into a genome browser like IGV.

To run:
Rscript GenerateTranslatome.R --TranslationCalls=path/to/translation_calls --NullDistribution=path/to/null_distribution --CandidateORFs=path/to/candidate_orfs

Options:
- --Output=path/to/output_folder: Designate output directory.
- --Threads=1: Specify thread count.
- --ExcludeChr=none: Chromosomes/contigs to exclude. Default is none. Example: --ExcludeChr=chr1,chr8,chrM
- --ExcludeOverlapGene=True: Exclude nORFs overlapping canonical genes on the same strand.
- --FDR=0.05: Define desired false discovery rate.
- --Scrambles=100: Set number of scrambles to calculate FDR. WIll run faster with less scrambles, if enough data is processed to maintain robustness.

------------------------------------------------------------------------------
