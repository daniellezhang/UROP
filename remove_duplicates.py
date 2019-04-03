#! /usr/bin/env python3

"""
Created on 12/02/19
Written by Danielle Zhang
Remove duplicate headers in SAM file
"""

import sys
import csv
import pysam



def remove_duplicate(sam_file):
    f = open(sam_file,"r")
    reader = csv.reader(f)
    header_lst = []
    names = []
    lengths = []
    header = {}

    #create a list of unique headers
    for row in reader:
        if not row or row[0][0:1] != '@':
            break
        elif row[0] not in header_lst:
            header_lst.append(row[0])

    f.close()


    #create a multi-level dictionary for header
    #extract name and length from header
    for h in header_lst:
        type = h.split('\t')[0][1:]
        if header.get(type,-1) == -1:
            header[type] = []
        line = {}
        line[h.split('\t')[1][:2]] = h.split('\t')[1][3:]
        length = h.split('\t')[2].split(':')
        line[h.split('\t')[2].split(':')[0]] = int(h.split('\t')[2].split(':')[1])
        if line not in header[type]:
            header[type].append(line)
            print(line)

    #create outfile with no duplicate header
    outfile = pysam.AlignmentFile("duplicates_removed.sam","w", header = header)
    samfile = pysam.AlignmentFile(sam_file,"r")
    #add alignments to outfile
    for read in samfile.fetch():
        outfile.write(read)

    samfile.close()
    outfile.close()

def main():
    #check if the arguments are valid
    assert len(sys.argv) == 2, "need to include a SAM file to be open"
    assert sys.argv[1].split('.')[-1]=='sam', "Input file need to be SAM format"
    sam_file = sys.argv[1]
    remove_duplicate(sam_file)


if __name__ == "__main__":
    main()
