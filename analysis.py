#! usr/bin/env python3

"""created on 19/03/19
written by Danielle Zhang
call other scripts to generate analysis for the alignment with the given\
fastq and SAM file"""

import sys
from count_lines import count
from uniquely_mapped_read import uniquely_mapped

def main():
    #check if the arguments are valid
    assert len(sys.argv) == 3, "Need to include the fastq file and SAM file\
    to be analysed"
    fastq_file = sys.argv[1]
    sam_file = sys.argv[2]

    count(fastq_file)
    uniquely_mapped(sam_file)


if __name__ == "__main__":
    main()
