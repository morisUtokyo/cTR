## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats from the same locus. It outputs a set of representative reads that are highly similar to each other. 

## Usage

cTR [-acs] [-f fasta file] [-d output directory] <fasta file name>
  
* -f: Input a fasta file name or a list of loci
* -d: Output the results to the specified directory 
* -a: Output a detailed analysis of groups
* -c: Output a fasta file of representative reads of groups

## Outline of the algorithm  
* We divide a set of reads from one individual into one or two alleles (haplotypes). For this, we measure the distance of any pair of reads, and group them into one or two groups such that the distance of any pair of elements is smaller than a given threshold. Assuming that the error rate of HiFi reads ranges from 0.1% to 0.5%, we set the maximum threshold to 1% (#define MAX_DIFF_RATIO  0.01 in clusterTR.h), which allows us to identify one or two groups in typical cases except for reads with many sequencing errors. From each group, we select a representative read that is the centroid in the group that minimizes the sum of the distances from it to the other reads. 

* We then collect a set of representative reads from individual genomes at each locus that are likely to have identical haplotypes. We then cluster the set by measuring their distances and grouping them into clusters.
  
* For calculating the distance beween a pair of DNA sequences, we align them using global alignment of KSW2, a library that aligns a pair of reads very efficiently and implements a suit of state-of-the-art dynamic programming techniques. For details, see:

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100

* For grouping reads, we used the Saitou-Nei neighbor-joining algorithm: See 
> * N Saitou, M Nei (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees., *Molecular Biology and Evolution*, Volume 4, Issue 4, Pages 406–425
