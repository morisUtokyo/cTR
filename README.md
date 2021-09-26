## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats from the same locus. It outputs a set of non-redundant representative reads that are dis-similar each other. 

## Usage

cTR [-ac] [-f input fasta file] [-d output directory]
  
* -f: Input a fasta file name
* -a: Output a detailed analysis of groups
* -c: Output a fasta file of representative reads of groups
* -d: Output the results (detailed analysis of groups and fasta files of representative reads) to the specified directory

For example, 

> cTR -ac -f 12345.fasta -d clustering_result

outputs

> 12345_analysis.txt and 12345_rep.fasta

to the directory named clustering_result.

## Outline of the algorithm  
* First, this program divides a set of PacBio HiFi reads collected from one individual into one or two alleles (haplotypes). For this, we measure the distance of any pair of reads, and group them into one or two groups such that the distance of any pair of elements is smaller than a given threshold. Assuming that the error rate of PacBio HiFi reads ranges from 0.1% to 0.5%, we set the maximum threshold to 1% (#define MAX_DIFF_RATIO  0.01 in clusterTR.h), which allows us to identify one or two groups typically except for reads with many sequencing errors. From each group, we select a representative read that is the centroid in the group that minimizes the sum of the distances from it to the other reads. Overall, one or two representative reads are output from one individual at each locus. 

* Second, this program clusters a set of representative reads from a number of individual genomes at each locus into groups such that the centroid of each group represents one tandem repeat haplotype.
  
* For calculating the distance beween a pair of DNA sequences, we align them using global alignment of KSW2, a library that aligns a pair of reads very efficiently and implements a suit of state-of-the-art dynamic programming techniques. For details, see:

> Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
>
> Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100

* For grouping reads, we used the Saitou-Nei neighbor-joining algorithm: See 
> N Saitou, M Nei (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees., *Molecular Biology and Evolution*, Volume 4, Issue 4, Pages 406–425
