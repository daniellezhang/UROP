#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 11:06:29 2019
Written by Danielle Zhang
test whether every sequence's quality string length equals to the actual
sequence length
"""
import sys
from Bio import SeqIO



#check if the arguments are valid
assert len(sys.argv) == 2, "need to include file to be open"

openfile = str(sys.argv[1])
iterator = SeqIO.parse(openfile, "fastq")

for read in iterator:
    """print("seq length= %d,quality_length=%d \n"%(len(read.seq),len(read.letter_annotations)))"""
    if len(read.seq) != len(read.letter_annotations['phred_quality']):
       print(read.id, end = '\n')
    


