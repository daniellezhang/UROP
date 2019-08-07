#!/usr/bin/env python3

"""
created on 18/07/19 by Danielle Zhang
extract a file of alignment information from the given BAM file"""

import pysam
import sys

def extract(filename):
    alignments = pysam.AlignmentFile(filename, 'rb')
    outfilename = filename.split('.')[0]+"alignment_info"
    outfile = open(outfilename,'w')

    for read in alignments.fetch(until_eof = True):
        outfile.write("%s,%s,%s,%s\n"%(read.query_name.split('/')[0],
        read.reference_name,read.reference_start,read.cigarstring))
    alignments.close()
    outfile.close()

def main():
    assert len(sys.argv) == 2, "need to include a BAM file"

    filename = sys.argv[1]
    extract(filename)

if __name__ == "__main__":
    main()
