## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats from the same locus. It outputs a set of representative reads that are highly similar to each other. 



cTR use KSW2, a library to align a pair of reads that implements dynamic programming. For details, see:

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100.

## Usage

cTR [-acs] [-f fasta file] [-d output directory] <fasta file name>
-f: Input a fasta file name or a list of loci
-d: Output the results to the specified directory 
-a: Output a detailed analysis of groups
-c: Output a fasta file of representative reads of groups
