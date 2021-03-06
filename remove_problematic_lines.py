#!/usr/bin/env python3

'''created on 05/02/19
Written by Danielle Zhang
script to remove the given problematic sequences and their neighbours'''

import sys
from subprocess import call


def remove(prob_seq, filename):

    n_neighbour_lines = 8
    line_list = []
    total_line = 0

    f = open(filename,'r')

    header = f.readline()
    seq = f.readline()
    plus = f.readline()
    quality_string = f.readline()

    #find the location of the problematic seqs
    while header and seq and plus and quality_string:
        seq_id = header.split(' ')[0][1:]
        if seq_id in prob_seq:
            line_list.append(total_line + 1)
        total_line += 4
        header = f.readline()
        seq = f.readline()
        plus = f.readline()
        quality_string = f.readline()

    f.close()
    line_list.sort()

    #create new fastq file with problematic seqs and their neighbours removed
    command = "awk '"
    current_line = 1
    for i in line_list:
        end = i - n_neighbour_lines-1
        if end > current_line:
            command += format("NR==%d,NR==%d;"%(current_line,end))
        current_line = i + n_neighbour_lines

    if current_line < total_line:
        command += format("NR==%d,NR==%d;"%(current_line,total_line))
    command +=format("' %s >./problematic_line_removed.fastq"%filename)

    call(command,shell = True)
    #check(prob_seq)

#function to check if all problematic lines have been removed
def check(prob_seq):
    f = open("./problematic_line_removed.fastq",'r')
    header = f.readline()
    seq = f.readline()
    plus = f.readline()
    quality_string = f.readline()
    while header and seq and plus and quality_string:
        seq_id = header.split(' ')[0][1:]
        if seq_id in prob_seq:
            print("wrong")
        header = f.readline()
        seq = f.readline()
        plus = f.readline()
        quality_string = f.readline()
    f.close()

def main():
    #check if the arguments are valid
    assert len(sys.argv) >= 3, "need to include file to be open and the sequence id \
    to be checked"

    #read lists of problematic seq id from the command line
    prob_seq = sys.argv[2:]
    filename = sys.argv[1]
    remove(prob_seq, filename)

if __name__ == "__main__":
    main()
