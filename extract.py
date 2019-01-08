#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Written by Danielle Zhang on 8/1/19
extracting sequences from reference genome and create a new fastq file
"""

import sys
from Bio import SeqIO

#check if the arguments are valid
assert len(sys.argv) < 3, "need to include a reference genome file \
and at least one transcript ID"

#list of transcript IDs
id = sys.argv[2:len(sys.argv)]

openfile = str(sys.argv[1])
iterator = SeqIO.parse(openfile, "fasta")
output = []

for record in iterator:
    if record.id in id:
        output.append(record)


SeqIO.write(output,"extract.fastq","fastq")