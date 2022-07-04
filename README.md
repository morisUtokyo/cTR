## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats from the same locus. It outputs a set of representative reads that are highly similar to each other. 

cTR use KSW2, a library to align a pair of reads that implements dynamic programming. Obtain a copy of the KSW2 program from:

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
