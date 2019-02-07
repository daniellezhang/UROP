#!/bin/bash

#created on 06/02/19
#written by Danielle Zhang
#Bash script to extract primary and secondary alignments respectively from given SAM file

#check if argument number is valid
if [[ $# -ne 1 ]]; then
	echo "Need to include a SAM file"
	exit
fi

filename=$1

#extract alignments from the input SAM file
#aligned reads with no supplementary alignment
{ samtools view -H $filename ; samtools view -G 4 $filename; } > temp.sam
{ samtools view -H $filename ; samtools view -G 2048 temp.sam; } > temp2.sam
rm -f temp.sam
#secondary alignment
{ samtools view -H $filename ; samtools view -f 256 temp2.sam; } > secondary.sam
#primary alignment
{ samtools view -H $filename ; samtools view -G 256 temp2.sam; } > primary.sam
rm -f temp2.sam

#extract read id from primary.sam and secondary.sam
samtools view primary.sam | cut -f 1 > primary.txt
samtools view secondary.sam | cut -f 1 > secondary.txt
