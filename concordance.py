#! /usr/bin/env python3

"""created on 31/07/19 by Danielle Zhang
concordance assessment of two alignment files
input files needed to be in BAM formats and are sorted by read names
"""


"""process CIGAR string to a list of aligned bases and the corresponding
position in reference sequence"""
def process_CIGAR(CIGAR_str):
    prev = 0
    current_pos = 0
    aligned_lst = []
    for i in range(len(CIGAR_str)):
        if CIGAR_str[i].isalpha():
            n = int(CIGAR_str[prev:i])
            if CIGAR_str[i] == 'M':
                aligned_lst.append((current_pos,current_pos+n))
            current_pos += n
            prev = i+1

    return aligned_lst


"""take two aligned bases lists to calculate concordance stats"""
def concordance_stat(aligned_lst1, aligned_lst2, ref_start1, ref_start2):
    aligned_by_both = 0
    uniquely_aligned_1 = 0
    uniquely_aligned_2 = 0
    unaligned_1 = 0
    unaligned_2 = 0
    unaligned_by_both = 0
    current_seq_pos1 = 0
    current_seq_pos2 = 0
    ref_pos = ref_start1


    #calculate the number of uniquely aligned bases when the leftmost coordinates are different
    if ref_start1 < ref_start2:
        i = 0
        while ref_pos<ref_start2:
             aligned_segment = aligned_lst1[i]
             aligned_bases = aligned_segment[1]-aligned_segment[0]
             if aligned_bases +ref_pos < ref_start2:
                 uniquely_aligned_1 += aligned_bases
                 ref_pos = ref_start1 + aligned_bases
                 i+=1
             else:
                 uniquely_aligned_1 += ref_start2 - ref_pos
                 aligned_lst1[i] = (ref_start2,aligned_lst1[i][1])
                 ref_pos = ref_start2


    #calculate concordnace stats from the same coordinates in the reference
    j = 0
    while i < len(aligned_lst1) and j < len(aligned_lst2):



    print(uniquely_aligned_1, uniquely_aligned_2)



def main():
    sample1 = "11S3M1D1I15M"
    sample2 = "5S7M3D10M"
    aligned_lst1 = process_CIGAR(sample1)
    aligned_lst2 = process_CIGAR(sample2)
    #ensure that the first list of aligned info has the leftest leftmost coordinates
    concordance_stat(aligned_lst1, aligned_lst2, 3,7)

if __name__ == "__main__":
    main()
