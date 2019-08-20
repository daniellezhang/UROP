#! /usr/bin/env python3

"""created on 31/07/19 by Danielle Zhang
concordance assessment of two alignment files
input files needed to be in BAM formats and are sorted by read names
"""
import sys


"a list of op that consumes query sequence"
consume_query_lst = ['M','I','S','=','X']

"a list of op that consume reference"
consume_ref_lst = ['D','N','M','=','X']

"a class to store the information extracted from CIGAR string"
class CIGAR_info(object):
    def __init__(self, query_lst, ref_lst):
        self.query_lst = query_lst
        self.ref_lst = ref_lst

"""process CIGAR string to a list of aligned bases and the corresponding
position in reference sequence"""
def process_CIGAR(CIGAR_str, ref_start):
    prev = 0
    current_query = 1
    current_ref = ref_start
    query_lst = []
    ref_lst = []
    #read CIGAR string by character
    for i in range(len(CIGAR_str)):
        #find a CIGAR operation
        if CIGAR_str[i].isalpha() or CIGAR_str[i] == '=':
            # number of bases for this CIGAR operation
            n = int(CIGAR_str[prev:i])
            op = CIGAR_str[i]
            #a match. add the segment of query and segment of ref to the lists
            if op == 'M':
                query_lst.append((current_query, current_query + n))
                ref_lst.append((current_ref, current_ref + n))
            #update the position in query sequence
            if op in consume_query_lst:
                current_query += n
            #update the position in reference sequence
            if op in consume_ref_lst:
                current_ref += n
            prev = i+1

    return CIGAR_info(query_lst,ref_lst)


"""take two aligned bases lists to calculate concordance stats
assume the reads mapped to the same reference segments"""
def concordance_stat(aligned1, aligned2, query_length):
    aligned_by_both = 0
    uniquely_aligned_1 = 0
    uniquely_aligned_2 = 0
    unaligned_1 = 0
    unaligned_2 = 0
    seq_start = min(aligned1.query_lst[0][0],aligned2.query_lst[0][0])
    seq_end = max(aligned1.query_lst[-1][1],aligned2.query_lst[-1][1])
    current_seq1 = aligned1.query_lst[0]
    current_seq2 = aligned2.query_lst[0]
    current_ref1 = aligned1.ref_lst[0]
    current_ref2 = aligned2.ref_lst[0]
    i = seq_start
    j = 0
    k = 0
    unaligned_by_both = seq_start + query_length - seq_end

    #check the overlap aligned bases base by base
    while i < seq_end and j < len(aligned1.query_lst) and k < len(aligned2.query_lst):
        #print(i,j,k)
        current_seq1 = aligned1.query_lst[j]
        current_seq2 = aligned2.query_lst[k]
        current_ref1 = aligned1.ref_lst[j]
        current_ref2 = aligned2.ref_lst[k]
        is_mapped1 = 0
        is_mapped2 = 0
        #base is mapped by both aligners
        if i == current_seq1[0] and i == current_seq2[0]:
            is_mapped1 = 1
            is_mapped2 = 1
            #base is mapped to the same position by both aligners
            if current_ref1[0] == current_ref2[0]:
                aligned_by_both += 1
            #base is mapped to different position by both aligners
            else:
                uniquely_aligned_1 += 1
                uniquely_aligned_2 += 1
        #base is only mapped by aligner 1
        elif i == current_seq1[0]:
            uniquely_aligned_1 += 1
            is_mapped1 = 1
        #base is only mapped by aligner 2
        elif i == current_seq2[0]:
            uniquely_aligned_2 += 1
            is_mapped2 = 1

        #update positions in query and reference
        if is_mapped1 == 1:
            if i+1 == current_seq1[1]:
                j += 1
            else:
                aligned1.query_lst[j] = (i+1, current_seq1[1])
                aligned1.ref_lst[j] = (current_ref1[0]+1, current_ref1[1])
        if is_mapped2 == 1:
            if i+1 == current_seq2[1]:
                k += 1
            else:
                aligned2.query_lst[k] = (i+1, current_seq2[1])
                aligned2.ref_lst[k] = (current_ref2[0]+1, current_ref2[1])

        if is_mapped1 == 0 and is_mapped2 == 0:
            unaligned_by_both += 1
        elif is_mapped1 == 0:
            unaligned_1 += 1
        elif is_mapped2 == 0:
            unaligned_2 += 1

        i += 1

    #one alignment has ended. continue calculate concordance stat
    #by base for the remainig alignment
    if j >= len(aligned1.query_lst):
        is_mapped1 = 0
        while  i < seq_end:
            #print(i,j,k)
            current_seq2 = aligned2.query_lst[k]
            current_ref2 = aligned2.ref_lst[k]
            is_mapped2 = 0
            if i == current_seq2[0]:
                is_mapped2 = 1
                uniquely_aligned_2 += 1

            if is_mapped2 == 1:
                unaligned_1 += 1
                if i+1 == current_seq2[1]:
                    k += 1
                else:
                    aligned2.query_lst[k] = (i+1, current_seq2[1])
                    aligned2.ref_lst[k] = (current_ref2[0]+1, current_ref2[1])
            else:
                unaligned_by_both += 1
            i += 1
    else:
        is_mapped2 = 0
        while  i < seq_end:
            #print(i,j,k)
            current_seq1 = aligned1.query_lst[j]
            current_ref1 = aligned1.ref_lst[j]
            is_mapped1 = 0
            if i == current_seq1[0]:
                is_mapped1 = 1
                uniquely_aligned_1 += 1

            if is_mapped1 == 1:
                unaligned_2 += 1
                if i+1 == current_seq1[1]:
                    j += 1
                else:
                    aligned1.query_lst[j] = (i+1, current_seq1[1])
                    aligned1.ref_lst[j] = (current_ref1[0]+1, current_ref1[1])
            else:
                unaligned_by_both += 1
            i += 1


    '''print("unaligned_by_both:%d, unaligned_1: %d, unaligned_2:%d,\
    uniquely_aligned_1: %d, uniquely_aligned_2: %d, aligned_by_both:%d"\
    %(unaligned_by_both,unaligned_1,unaligned_2,uniquely_aligned_1,\
    uniquely_aligned_2,aligned_by_both))'''

    return (unaligned_by_both,unaligned_1,unaligned_2,\
    uniquely_aligned_1,uniquely_aligned_2,aligned_by_both)

'''take two aligned base lists to calculate concordance stat
assume that the the reads mapped to different reference segments'''
def concordance_stat2(aligned1, aligned2, query_length):
    aligned_by_both = 0
    uniquely_aligned_1 = 0
    uniquely_aligned_2 = 0
    unaligned_1 = 0
    unaligned_2 = 0
    seq_start = min(aligned1.query_lst[0][0],aligned2.query_lst[0][0])
    seq_end = max(aligned1.query_lst[-1][1],aligned2.query_lst[-1][1])
    current_seq1 = aligned1.query_lst[0]
    current_seq2 = aligned2.query_lst[0]
    i = seq_start
    j = 0
    k = 0
    unaligned_by_both = seq_start + query_length - seq_end


    #check the overlap aligned bases base by base
    while i < seq_end and j < len(aligned1.query_lst) and k < len(aligned2.query_lst):
        #print(i,j,k)
        current_seq1 = aligned1.query_lst[j]
        current_seq2 = aligned2.query_lst[k]
        is_mapped1 = 0
        is_mapped2 = 0

        #base is mapped by aligner 1
        if i == current_seq1[0]:
            uniquely_aligned_1 += 1
            is_mapped1 = 1
        #base is mapped by aligner 2
        if i == current_seq2[0]:
            uniquely_aligned_2 += 1
            is_mapped2 = 1

        #update positions in query
        if is_mapped1 == 1:
            if i+1 == current_seq1[1]:
                j += 1
            else:
                aligned1.query_lst[j] = (i+1, current_seq1[1])
        if is_mapped2 == 1:
            if i+1 == current_seq2[1]:
                k += 1
            else:
                aligned2.query_lst[k] = (i+1, current_seq2[1])

        if is_mapped1 == 0 and is_mapped2 == 0:
            unaligned_by_both += 1
        elif is_mapped1 == 0:
            unaligned_1 += 1
        elif is_mapped2 == 0:
            unaligned_2 += 1

        i += 1

    #one alignment has ended. continue calculate concordance stat
    #by base for the remainig alignment
    if j >= len(aligned1.query_lst):
        is_mapped1 = 0
        while  i < seq_end:
            current_seq2 = aligned2.query_lst[k]
            is_mapped2 = 0
            if i == current_seq2[0]:
                is_mapped2 = 1
                uniquely_aligned_2 += 1

            if is_mapped2 == 1:
                unaligned_1 += 1
                if i+1 == current_seq2[1]:
                    k += 1
                else:
                    aligned2.query_lst[k] = (i+1, current_seq2[1])
            else:
                unaligned_by_both += 1
            i += 1
    else:
        is_mapped2 = 0
        while  i < seq_end:
            current_seq1 = aligned1.query_lst[j]
            is_mapped1 = 0
            if i == current_seq1[0]:
                is_mapped1 = 1
                uniquely_aligned_1 += 1

            if is_mapped1 == 1:
                unaligned_2 += 1
                if i+1 == current_seq1[1]:
                    j += 1
                else:
                    aligned1.query_lst[j] = (i+1, current_seq1[1])
            else:
                unaligned_by_both += 1
            i += 1

    return (unaligned_by_both,unaligned_1,unaligned_2,\
    uniquely_aligned_1,uniquely_aligned_2,aligned_by_both)

#driver function for running the concordance analysis
def concordance_driver():
    assert len(sys.argv) == 2, "Need to include an input file"
    mapping_info = sys.argv[1]

    '''read the mapping information from the given csv file
    headers: row, qname, qlength, ref1, ref_start1, CIGAR1, ref2, ref_start2, CIGAR2 '''
    f1 = open(mapping_info, 'r')
    f1.readline()


    '''output a csv file of concordance stat, one query per line'''
    f2 = open("concordance_stat.csv", 'w+')
    f2.write("qname, same_ref, unaligned_both, unaligned_1, unaligned_2, aligned_1,aligned_2, aligned_both\n")
    for query in f1.readline():
        query = f1.readline()[:-1].split(',')
        qname = query[1][1:-1]
        qlength = int(query[2])
        ref1 = query[3][1:-1]
        ref_start1 = int(query[4])
        CIGAR1 = query[5][1:-1]
        ref2 = query[6][1:-1]
        ref_start2 = int(query[7])
        CIGAR2 = query[8][1:-1]
        aligned1 = process_CIGAR(CIGAR1, ref_start1)
        aligned2 = process_CIGAR(CIGAR2, ref_start2)
        same_ref = 0
        if ref1 == ref2:
            stat = concordance_stat(aligned1, aligned2, qlength)
            same_ref = 1
        else:
            stat = concordance_stat2(aligned1, aligned2, qlength)
        str_line = str(stat)[1:-1]
        line = format("%s,%d,%s\n"%(qname, same_ref, str_line))
        f2.write(line)


    f1.close()
    f2.close()

def main():
    #test
    '''sample1 = "30M1D3M1D5M1I2M1I6M1I4M1077N1M1D56M1I46M3D9M415N34M2D1M1I32M2D24M2D4M1D10M2D90M2D28M1D26M2D9M1D13M1D40M1D44M5D11M2I2M1D2M2D22M94S"
    sample2 = "1S28M1D4M1D5M1I2M1I6M1I4M1077N1M1D55M1I48M2D1M1D7M415N34M2D1M1I34M1D1M1D25M3D11M2D46M1I3M1D23M1D4M1I12M2D29M1D24M2D13M1D13M1D38M1D48M5D5M2I2M1D5M1D3M1D16M5D3M2D5M1D4M1D3M1D5M3D1M73S"
    aligned1 = process_CIGAR(sample1, 12975116)
    aligned2 = process_CIGAR(sample2, 12975117)
    concordance_stat(aligned2, aligned1,655)
    sample1 = "5S3M1D1I15M"
    sample2 = "5S7M3D10M"
    aligned1 = process_CIGAR(sample1, 5)
    print(aligned1.ref_lst, aligned1.query_lst)
    aligned2 = process_CIGAR(sample2, 5)
    print(aligned2.ref_lst, aligned2.query_lst)
    #ensure that the first aligned info has the leftest leftmost coordinates in ref
    concordance_stat(aligned1, aligned2,40)'''

    concordance_driver()


if __name__ == "__main__":
    main()
