## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats from the same locus. It feeds a fasta file of mosaic tandem repeats at a focal locus and assumes that each read is annotated with the individual identifier and the read identifier. For example, in the following fasta file, the first annotation "0,0,(ACC)32(GTTT)30" shows the individual ID is 0, the read ID is 0, and the pattern of the mosaic tandem repeat is (ACC)32(GTTT)30, which means ACC occurs 32 times and GTTT 30 times. The string in the second row is identical to the pattern.

> 0,0,(ACC)32(GTTT)30
ACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTT
> 0,1,(ACC)32(GTTT)30
ACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCACCGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTT

cTR outputs a set of representative reads that are highly similar to each other. 

cTR uses KSW2, a library to align a pair of reads that implements dynamic programming. Obtain a copy of the KSW2 program from:

https://github.com/lh3/ksw2

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100.

For your convenience, a copy of the KSW2 program is placed on this github.

## Usage

cTR [-ac] [-f fasta file] [-d output directory]
* -f: Input a fasta file name without .fasta
* -i: Input directory
* -d: Output the results to the specified directory 
* -a: Output a detailed analysis of groups
* -c: Output a fasta file of representative reads of groups

## Test

To check how cTR works, go to the directory test and execute test.sh, which automatically generates a sample fasta file and applies cTR to the file. In test.sh, one can change the values of a variety of parameters for creating various samples. Initial values are set to these parameters, and test.sh generates the following files:

* ACC_GTTT.fasta: A sample fasta file with a total of 400 reads (two reads for two haplotypes of 100 individuals). Mosaic tandem repeats of the form (ACC)m(GTTT)n, where m and n range from 30 to 33.
* ACC_GTTT_rep.fasta: A fasta file after application of cTR to the sample fasta file. It contains the representative (centroid) mosaic tandem repeat in each group, and the annotation of the representative shows the number of elements in each cluster (e.g., GroupSize = 46), the diameter and radius of the cluster (e.g., Diameter = 3), the individual and read IDs of the centroid with its pattern (e.g., CentroidReadName =  2,2,(ACC)31(GTTT)32), and the length of the centroid repeat (e.g., CentroidReadLength = 221).
* ACC_GTTT_analysis.fasta: Describes the details of each centroid and its members, and the neighbor joining (NJ) tree of clusters.
* ACC_GTTT_table.fasta: The first column shows each mosaic tandem repeat with its individual and read IDs, and its repeat pattern. The second column display the cetroid of the group that the repeat in the first column belongs to.
