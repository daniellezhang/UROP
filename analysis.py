#! usr/bin/env python3

"""created on 19/03/19
written by Danielle Zhang
call other scripts to generate analysis for the alignment with the given\
fastq and SAM file"""

import sys
from count_lines import count
from uniquely_mapped_read import uniquely_mapped
from subprocess import call

def main():
    #check if the arguments are valid
    assert len(sys.argv) ==2, "Need to include the SAM file to be analysed"
    sam_file = sys.argv[1]

    uniquely_mapped(sam_file)
    print("number of supplementary alignment:")
    call("samtools view -c -f 2048 "+sam_file, shell=True)
    print("number of unaligned reads:")
    call("samtools view -c -f 4 "+sam_file, shell=True)
    print("number of secondary alignment:")
    call("samtools view -c "+"./secondary.sam", shell=True)
    print("number of primary alignment:")
    call("samtools view -c "+"./primary.sam", shell=True)
    call("bash ./aligned_base_count.sh " + "./primary.sam", shell = True)
    call("bash ./remove_supplementary_files.sh",shell=True)


if __name__ == "__main__":
    main()
