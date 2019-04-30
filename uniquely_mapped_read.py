#! usr/bin/env python3

"""created on 07/02/19
written by Danielle Zhang
calculate the number of unique reads
given the SAM file"""

import sys
from subprocess import call

def uniquely_mapped(sam_file):
    #extract primary and secondary alignments' read ids
    call("./SAM_file_manipulation.sh " + sam_file, shell = True)


    try:
        #build a dictionary of primary alignments' ids and count
        #the number of their occurence
        id_dict = {}
        f = open('./primary.txt','r')
        for id in f.readlines():
            if not id_dict.get(id):
                id_dict[id] = 1
            else:
                id_dict[id] += 1
        f.close()

        #count the number of occurence of secondary alignments' id
        f = open('./secondary.txt','r')
        for id in f.readlines():
            if not id_dict.get(id):
                id_dict[id] = 1
            else:
                id_dict[id] += 1
        f.close()

        #count the number of uniquely mapped read
        n_unique_mapped = 0
        for key in id_dict.keys():
            if id_dict[key] == 1:
                n_unique_mapped += 1

        print("Number of uniquely mapped read: %d\n"%(n_unique_mapped))


    except FileNotFoundError:
       print("primary and secondary alignments' read id files not found")
       sys.exit()

def main():
        #check if the arguments are valid
        assert len(sys.argv) == 2, "need to include a SAM file to be open"

        sam_file = sys.argv[1]
        uniquely_mapped(sam_file)



if __name__ == "__main__":
    main()
