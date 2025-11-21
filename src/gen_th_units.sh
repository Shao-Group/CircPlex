#!/bin/bash

tool_dir=/data/tkz5115/repeat_project/tools/TideHunter-v1.5.5/bin
input_file=$1
th_file=$2
output_file=$3
estm_length_file="/data/tkz5115/rev_comp_reads_analysis/unit_length_estimation_bucket/SRR10612066_complex_k10"

echo $th_file

$tool_dir/TideHunter -f 2 $input_file > $th_file

# Compile the program
g++ -o th_to_fasta th_to_fasta.cpp -std=c++11

if [ $? -eq 0 ]; then
    echo "th_to_fasta compilation successful."

    echo $output_file
    
    ./th_to_fasta $input_file $th_file $estm_length_file $output_file
else
    echo "Compilation failed."
fi