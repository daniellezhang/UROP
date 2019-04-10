#!/bin/bash

#created on 06/02/19
#written by Danielle Zhang
#Bash script to count the number of aligned base in the given SAM/BAM file

#check if argument number is valid
if [[ $# -ne 1 ]]; then
	echo "Need to include a SAM file"
	exit
fi

filename=$1

#convert SAM file to BAM
samtools view -Sb $filename | samtools sort > sorted_aln.bam

#use bedtools to output alignment base information
bedtools genomecov -ibam sorted_aln.bam -bga > genomcov_output.bam
#count the number of aligned base
awk '{s+=$4}END{print "mapped read bases: " s}' genomcov_output.bam
