#! /usr/bin/env python3

"""
Created on 06/02/19
Written by Danielle Zhang
Change the header of the input fastq file to pacbio's format
"""

import sys

def rename(filename):
    newfile = open("renamed.fastq","wt")
    f = open(filename,'r')
    total_seq = 0

    header=f.readline()

    #parse the input to generate new file
    while header:
        total_seq += 1
        seq = f.readline()
        plus = f.readline()
        quality_string = f.readline()
        old_header = header.split(' ',1)
        #create new header
        new_header = format("%s/%d/0_%d %s"%(old_header[0],total_seq,len(seq)-1,old_header[1]))
        newfile.write(new_header)
        newfile.write(seq)
        newfile.write(plus)
        newfile.write(quality_string)
        header=f.readline()


    f.close()
    newfile.close()


def main():
    #check if the arguments are valid
    assert len(sys.argv) == 2, "Need to include fastq file to be renamed in the arguments"
    assert sys.argv[1].split('.')[-1]=='fastq', "Input file need to be fastq format"

    filename = sys.arv[1]
    rename(filename)


if __name__ == "__main__":
    main()
