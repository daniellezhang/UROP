#! /usr/bin/env python3

"""
Created on 12/02/19
Written by Danielle Zhang
Remove duplicate headers in SAM file
"""

import sys
import pysam
import csv

#check if the arguments are valid
assert len(sys.argv) == 2, "need to include a SAM file to be open"
assert sys.argv[1].split('.')[-1]=='sam', "Input file need to be SAM format"

f = open(sys.argv[1],"r")
reader = csv.reader(f)
header_lst = []
names = []
lengths = []

#create a list of unique headers
for row in reader:
    if not row or row[0][0:1] != '@':
        break
    elif row[0] not in header_lst:
        header_lst.append(row[0])

f.close()

#extract name and length from header
for h in header_lst:
    names.append(h.split('\t')[1][3:])
    lengths.append(int(h.split('\t')[2][3:]))

#create AlignmentHeader file
header = pysam.AlignmentHeader.from_references(reference_lengths = lengths,reference_names = names)
#create outfile with no duplicate header
outfile = pysam.AlignmentFile("duplicates_removed.sam","w", header = header)
samfile = pysam.AlignmentFile(sys.argv[1],"r")
#add alignments to outfile
for read in samfile.fetch():
    outfile.write(read)

samfile.close()
outfile.close()
