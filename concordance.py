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
    query_lst = []
    ref_lst = []
    ref_start = int(ref_start)
    current_ref = ref_start
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
def concordance_stat(aligned1, aligned2, query_length, same_ref):
    aligned_by_both = 0
    aligned_1= 0
    aligned_2 = 0
    aligned_different_pos = 0
    seq_start = min(aligned1.query_lst[0][0],aligned2.query_lst[0][0])
    seq_end = max(aligned1.query_lst[-1][1],aligned2.query_lst[-1][1])

    i = seq_start
    j = 0
    k = 0
    #unaligned bases at the start = seq_start-1
    #unaligned bases at the end = query_length-(seq_end-1)
    unaligned_by_both = seq_start + query_length - seq_end

    #check the overlap aligned bases base by base
    while i < seq_end:
        #print(i,j,k)
        is_mapped1 = 0
        is_mapped2 = 0

        if j < len(aligned1.query_lst):
            current_seq1 = aligned1.query_lst[j]
            current_ref1 = aligned1.ref_lst[j]
        elif j == len(aligned1.query_lst):
            current_seq1 = (-1,-1)

        if k < len(aligned2.query_lst):
            current_seq2 = aligned2.query_lst[k]
            current_ref2 = aligned2.ref_lst[k]
        elif k == len(aligned2.query_lst):
            current_seq2 = (-1,-1)
        #base is mapped by both aligners
        if i == current_seq1[0] and i == current_seq2[0]:
            is_mapped1 = 1
            is_mapped2 = 1
            #base is mapped to the same position by both aligners
            if current_ref1[0] == current_ref2[0] and same_ref == 1:
                aligned_by_both += 1
            #base is mapped to different position by both aligners
            else:
                aligned_different_pos += 1
        #base is only mapped by aligner 1
        elif i == current_seq1[0]:
            is_mapped1 = 1
            aligned_1 += 1
        #base is only mapped by aligner 2
        elif i == current_seq2[0]:
            is_mapped2 = 1
            aligned_2 += 1

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


        i += 1



    return (unaligned_by_both,aligned_by_both, aligned_different_pos, aligned_1, aligned_2)

#count the number of aligned bases
def aligned_base_count(CIGAR_str):
    aligned_count = 0
    prev = 0
    for i in range(len(CIGAR_str)):
        #find a CIGAR operation
        if CIGAR_str[i].isalpha() or CIGAR_str[i] == '=':
            # number of bases for this CIGAR operation
            n = int(CIGAR_str[prev:i])
            op = CIGAR_str[i]
            #a match. add the segment of query and segment of ref to the lists
            if op == 'M' or op == '=':
                aligned_count += n
            prev = i+1
    return aligned_count





#driver function for running the concordance analysis
def concordance_driver():
    assert len(sys.argv) == 2, "Need to include an input file"
    mapping_info = sys.argv[1]

    '''read the mapping information from the given csv file
    headers: row, qname, qlength, aligner1, ref1, ref_start1, CIGAR1, aligner2, ref2, ref_start2, CIGAR2 '''
    f1 = open(mapping_info, 'r')

    '''output a csv file of concordance stat, one query per line'''
    f2 = open("concordance_stat.csv", 'w+')
    f2.write("qname, qlength, unaligned_both, aligned_by_both, aligned_different_pos, aligned_1, aligned_2\n")
    for line in f1.readlines()[1:]:
        query = line[:-1].split(',')
        qname = query[1][1:-1]
        qlength = int(query[2])
        aligner1 = int(query[3])
        ref1 = query[4][1:-1]
        ref_start1 = query[5]
        CIGAR1 = query[6][1:-1]
        aligner2 = int(query[7])
        ref2 = query[8][1:-1]
        ref_start2 = query[9]
        CIGAR2 = query[10][1:-1]
        same_ref = 0
        #read is mapped by both aligners
        if aligner1 == 1 and aligner2 == 1:
            ref_start1 = int(ref_start1)
            ref_start2 = int(ref_start2)
            alignment1 = process_CIGAR(CIGAR1, ref_start1)
            alignment2 = process_CIGAR(CIGAR2, ref_start2)
            #read is mapped to the same reference
            if ref1 == ref2:
                same_ref = 1
            stat = concordance_stat(alignment1, alignment2, qlength, same_ref)
            str_line = str(stat)[1:-1]
            line = format("%s,%d,%s\n"%(qname, qlength, str_line))
        #read is mapped by only one aligner
        elif aligner1 == 1:
            aligned_base = aligned_base_count(CIGAR1)
            unaligned_both = qlength - aligned_base
            line = format("%s,%d,%d,%d,%d,%d,%d\n"%(qname, qlength, unaligned_both,0,0,aligned_base,0))
        elif aligner2 == 1:
            aligned_base = aligned_base_count(CIGAR2)
            unaligned_both = qlength - aligned_base
            line = format("%s,%d,%d,%d,%d,%d,%d\n"%(qname, qlength, unaligned_both,0,0,0,aligned_base))
        #read is unaligned by both aligners
        else:
            line = format("%s,%d,%d,0,0,0,0\n"%(qname,qlength, qlength))
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
