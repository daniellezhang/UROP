#1 /usr/bin/python
"""
Written by Danielle Zhang on 14/11/2018
create a shorter fastq file from the input fastq file for testing
"""

import sys
from Bio import SeqIO

#check if the arguments are valid
assert len(sys.argv) == 2, "need to include file to be open "

openfile = str(sys.argv[1])

#use the first 8 reads in the input file to create a new fastq file
n = 8
i = 0

input_seq_iterator = SeqIO.parse(openfile, "fastq")
output = []

for read in input_seq_iterator:
    i+=1
    if i == n:
        break
    output.append(read)


SeqIO.write(output, "test.fastq","fastq")
