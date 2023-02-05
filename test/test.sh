#!/bin/bash
# A bash script file for testing the program, cTR, that clusters reads in each fasta file
# Call bash test.sh

num_read_in_one_haplotype=2 # Number of reads in one of two haplotypes
num_individuals=100 # Number of individuals
num_units=2     # Number of different units in mosaic tandem repeats
unit1="ACCCC"     # String of the first unit
unit2="GTTTT"    # String of the second unit
# Mosaic tandem repeats of the above units whose frequency of occurrence is between the lower and upper limits below.
lower_bound_num_unit_occurrences=20
upper_bound_num_unit_occurrences=25

cTR=../cTR # executable modules
result=./  # the directory where results are output

sample=$unit1"_"$unit2

cd gendata; make clean; make; cd ..;

gendata/gen -h $num_read_in_one_haplotype -n $num_individuals -k $lower_bound_num_unit_occurrences -l $upper_bound_num_unit_occurrences -m $num_units $unit1 $unit2 > $sample.fasta

$cTR -ac -f $sample -d $result

exit 0
