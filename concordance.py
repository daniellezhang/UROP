#! /usr/bin/env python3

"""created on 31/07/19 by Danielle Zhang
concordance assessment of two alignment files
input files needed to be in BAM formats and are sorted by read names
"""

"a list of op that  consumes query sequence"
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
    current_query = 0
    current_ref = ref_start
    query_lst = []
    ref_lst = []
    for i in range(len(CIGAR_str)):
        if CIGAR_str[i].isalpha() or CIGAR_str[i] == '=':
            n = int(CIGAR_str[prev:i])
            op = CIGAR_str[i]
            #a match. add the segment of query and segment of ref to the lists
            if op == 'M':
                query_lst.append((current_query ,current_query + n))
                ref_lst.append((current_ref, current_ref + n))
            #update the position in query sequence
            if op in consume_query_lst:
                current_query += n
            #update the position in reference sequence
            if op in consume_ref_lst:
                current_ref += n
            prev = i+1
    print(query_lst, ref_lst)
    return CIGAR_info(query_lst,ref_lst)


"""take two aligned bases lists to calculate concordance stats
assume that the first list has the leftest leftmost """
def concordance_stat(aligned1, aligned2):
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
    unaligned_by_both = seq_start
    #check the overlap aligned bases base by base
    while i < seq_end and j < len(aligned1.query_lst) and k < len(aligned2.query_lst):
        current_seq1 = aligned1.query_lst[j]
        current_seq2 = aligned1.query_lst[k]
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
        elif i == current_seq1[0] and i != current_seq2[0]:
            uniquely_aligned_1 += 1
            is_mapped1 = 1
        #base is only mapped by aligner 2
        elif i != current_seq1[0] and i == current_seq2[0]:
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

    print("unaligned_by_both:%d, unaligned_1: %d, unaligned_2:%d,\
    uniquely_aligned_1: %d, uniquely_aligned_2: %d, aligned_by_both:%d"\
    %(unaligned_by_both,unaligned_1,unaligned_2,uniquely_aligned_1,uniquely_aligned_2,aligned_by_both))

    return


def main():
    sample1 = "5S3M1D1I15M"
    sample2 = "5S7M3D10M"
    aligned1 = process_CIGAR(sample1, 5)
    aligned2 = process_CIGAR(sample2, 7)
    #ensure that the first list of aligned info has the leftest leftmost coordinates in ref
    concordance_stat(aligned1, aligned2)

if __name__ == "__main__":
    main()
