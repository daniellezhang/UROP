#! /usr/bin/env python3

"""
Created on 12/02/19
Written by Danielle Zhang
Remove duplicate headers in SAM file
"""

import sys
from subprocess import call

#check if the arguments are valid
assert len(sys.argv) == 2, "need to include a SAM file to be open"
assert sys.argv[1].split('.')[-1]=='sam', "Input file need to be SAM format"

f = open(sys.argv[1],'r')
fout = open("duplicates_removed.sam", 'w')

header_tags = ["@HD", "@SQ", "@PG", "@RG", "@CO"]

headers = []
for line in f.readlines():
    #not a header, add to the output file
    if line.split(" ")[0] not in header_tags:
        fout.write(line)
    else:
        #add only unique headers to the output file
        if line not in headers:
            headers.append(line)
            fout.write(line)

f.close()
fout.close()
