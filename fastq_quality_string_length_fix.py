# -*- coding: utf-8 -*-

"""
Written by Danielle Zhang on 14/11/2018
fixing the quality string length in every header of fastq file 
to equal to the actual sequence length
"""

import sys
from Bio import SeqIO


#check if the arguments are valid
assert len(sys.argv) == 3, "need to include file to be open and the file\
to be created"

openfile = str(sys.argv[1])
writefile = str(sys.argv[1])
