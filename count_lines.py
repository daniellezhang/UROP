#! usr/bin/env python3

"""created on 19/03/19
written by Danielle Zhang
count the number of read in fastq file """

import sys

def count(filename):
    count = 0

    f = open(filename,'r')
    header=f.readline()

    #parse the input to generate new file
    while header:
        count += 1
        seq = f.readline()
        plus = f.readline()
        quality_string = f.readline()
        header=f.readline()

    f.close()
    print("number of read: %d"%(count))

def main():
    #check if the arguments are valid
    assert len(sys.argv) == 2, "Need to include fastq file to be renamed in the arguments"

    filename = sys.arv[1]
    count(filename)

if __name__ == "__main__":
    main()
