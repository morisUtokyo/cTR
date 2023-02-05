## Introduction

cTR is a tool for clustering PacBio HiFi reads with mosaic tandem repeats that are collected from multiple samples at the same locus. It outputs a set of representative reads that are highly similar to each other. 

cTR feeds such a fasta file that the annotation begins with a sample ID and a read ID. For example, the first row of ACCCC_GTTTT.fasta in the "test" directory shows that the sample and read IDs are 0 and 0, and the rest indicates mosaic tandem repeat pattern of the read string:

    > 0,0,(ACCCC)22(GTTTT)22
    ACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT
    > 0,1,(ACCCC)22(GTTTT)22
    ACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT

cTR clusters reads from an identical sample into one or two groups, and outputs the representative each group. If two haplotypes are homozygous, one group is output; if two haplotypes are heterozygous, two groups are output. In one group, the maximum Levenshtein distance between a pair of reads is 0.01 times the average length of the reads. Each of two groups also meets this condition. The default maximum value is set ot 0.01 because the maximum error rate of PacBio reads is less than 0.01 as of today. But it can be changed (say, 0.005) by redefining MAX_DIFF_RATIO in clusterTR.h.

Subsequently, cTR further clusters representative reads from all samples into groups so that the maximum Levenshtein distance between a pair of reads in a group is smaller than a threshold. The maximum threshold is set to 0.03 by default but can be also changed by redefining MAX_DIAMETER in cluster.h. cTR outputs a representative of each group with an annotation of the form:

    > GroupSize = N, Diameter = D, RadiusFromCentroid = R, CentroidReadName =  sampleID,readID, CentroidReadLength = L

For example, the first two rows of ACCCC_GTTTT_rep.fasta in the "test" directory is:

    > GroupSize = 39, Diameter = 10, RadiusFromCentroid = 10, CentroidReadName =  6,2,(ACCCC)21(GTTTT)23, CentroidReadLength = 220
    ACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCACCCCGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTTGTTTT

Each representative read of a sample is put into one of the groups in all samples. cTR outputs a table that describe this association. For example, the first five rows in ACCCC_GTTTT_hap_cent_table.txt are:

    0,0,(ACCCC)22(GTTTT)22     9,2,(ACCCC)22(GTTTT)23
    0,2,(ACCCC)24(GTTTT)20     25,2,(ACCCC)23(GTTTT)20
    1,0,(ACCCC)24(GTTTT)21     41,0,(ACCCC)24(GTTTT)22
    1,2,(ACCCC)24(GTTTT)20     25,2,(ACCCC)23(GTTTT)20
    2,0,(ACCCC)23(GTTTT)23     8,0,(ACCCC)24(GTTTT)23

cTR also generates a table that shows a list of representative reads in each group of all samples. For example, ACCCC_GTTTT_table.txt shows:

     6,2,(ACCCC)21(GTTTT)23     6,2,(ACCCC)21(GTTTT)23
     15,2,(ACCCC)21(GTTTT)23     6,2,(ACCCC)21(GTTTT)23
     17,0,(ACCCC)21(GTTTT)23     6,2,(ACCCC)21(GTTTT)23
     21,0,(ACCCC)21(GTTTT)23     6,2,(ACCCC)21(GTTTT)23
     26,2,(ACCCC)21(GTTTT)23     6,2,(ACCCC)21(GTTTT)23    

## Usage

cTR [-ac] [-f fasta file] [-d output directory]
* -f: Input a fasta file name without .fasta
* -i: Input directory
* -d: Output the results to the specified directory 
* -a: Output a detailed analysis of groups
* -c: Output a fasta file of representative reads of groups

For example, see test.sh in the test directory.


cTR use KSW2, a library to align a pair of reads that implements dynamic programming. Obtain a copy of the KSW2 program from:

https://github.com/lh3/ksw2

> * Suzuki, H. and Kasahara, M. (2018). Introducing difference recurrence relations for faster semi-global alignment of long sequences. *BMC Bioinformatics*, **19**:45.
> * Li, H (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, **34**:3094-3100.

For your convenience, a copy of the KSW2 program is placed on this github.




