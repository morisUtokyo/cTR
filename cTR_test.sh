#!/bin/bash
# A bash script file for clustering reads in each fasta file in the current directory


cTR=../v5/cTR           # executable modules
result=result_v5_0.01/  # the directory where results are output

echo "This program cluster tandem repeats (TRs) given in " $result

original_ifs=$IFS
IFS=$'\n'

file_names=($(eval ls $dir | grep 'fasta'))
echo "${file_names[@]}"

IFS=$original_ifs

for f in "${file_names[@]}"; do
    echo $f
    $cTR -ac -f $f -d $result
done

exit 0
