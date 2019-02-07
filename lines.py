#!/usr/bin/env python3

'''created on 05/02/19 
Written by Danielle Zhang
script to return the number of lines that the given read id is on '''

import sys
from Bio import SeqIO

#check if the arguments are valid
assert len(sys.argv) >= 3, "need to include file to be open and the sequence id \
to be checked"

#read lists of problematic seq id from the command line
prob_seq = sys.argv[2:]

total_line = 0
openfile = sys.argv[1]
seq_iterator = SeqIO.parse(openfile, "fastq")

for read in seq_iterator:
    if read.id in prob_seq:
       print("%s %d \n" %(read.id, total_line+1))
    total_line += 1
    
print(total_line)