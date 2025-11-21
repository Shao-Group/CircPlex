#!/bin/bash

input=$1
output=$2
minimap2_dir=$3
ref=/data/tkz5115/human_GRCh37_ensembl/Homo_sapiens.GRCh37.dna.primary_assembly.fa
anno=/data/tkz5115/human_GRCh37_ensembl/Homo_sapiens.GRCh37.87.bed #needs bed file
index_dir=/data/tkz5115/human_GRCh37_ensembl/minimap2-index

# # indexing
# mkdir -p $index_dir
# $minimap2_dir/minimap2 -d $index_dir/human_ref.mmi $ref

#mapping long ont reads
# $minimap2_dir/minimap2 --junc-bed $anno -t 20 -Y -ax splice -ub --MD --eqx $ref $input > $output
# $minimap2_dir/minimap2 --junc-bed $anno -t 30 -Y -ax splice -ub --MD $ref $input > $output
$minimap2_dir/minimap2 --junc-bed $anno -t 10 -Y -ax splice -ub --MD $index_dir/human_ref.mmi $input > $output