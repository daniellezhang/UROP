# -*- coding: utf-8 -*-

"""
Written by Danielle Zhang on 14/11/2018
fixing the quality string length in every header of fastq file 
to equal to the actual sequence length
"""

import sys
from Bio import SeqIO



#check if the arguments are valid
assert len(sys.argv) == 2, "need to include file to be open"

openfile = str(sys.argv[1])

input_seq_iterator = SeqIO.parse(openfile, "fastq")
output = []

for read in input_seq_iterator:
   output.append(read)
   description = output[len(output)-1].description
   i = description.find("ch=")
   j = description.find(' ', i)
   seq_len = len(read.seq)
   #replace the quality string length by the actual sequence length
   description = description[0:i]+"ch="+str(seq_len)+description[j:]

   output[len(output) -1].description = description


#create new fastq with the correct quality string length
SeqIO.write(output,"output.fastq", "fastq")