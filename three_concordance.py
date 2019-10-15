import sys
from concordance import process_CIGAR, concordance_stat,aligned_base_count



"""take three aligned bases lists to calculate concordance stats
assume the reads are mapped by three aligners"""
def concordance_stat3(aligned1, aligned2, aligned3, query_length, ref1, ref2, ref3):
    same_pos_3 = 0
    same_pos_1_2_diff_3  = 0
    same_pos_2_3_diff_1 = 0
    same_pos_1_3_diff_2 = 0
    same_pos_1_2_unmapped_3 = 0
    same_pos_2_3_unmapped_1 = 0
    same_pos_1_3_unmapped_2 = 0
    diff_pos_3 = 0
    diff_pos_1_2_unmapped_3 = 0
    diff_pos_2_3_unmapped_1 = 0
    diff_pos_1_3_unmapped_2 = 0
    mapped_only_1 = 0
    mapped_only_2 = 0
    mapped_only_3 = 0

    seq_start = min(aligned1.query_lst[0][0],aligned2.query_lst[0][0],aligned3.query_lst[0][0])
    seq_end = max(aligned1.query_lst[-1][1], aligned2.query_lst[-1][1],aligned3.query_lst[-1][1])

    i = seq_start
    j = 0
    k = 0
    l = 0
    all_unmapped = seq_start + query_length - seq_end


    #check the overlap aligned bases base by base
    while i < seq_end:
        is_mapped1 = 0
        is_mapped2 = 0
        is_mapped3 = 0

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

        if l < len(aligned3.query_lst):
            current_seq3 = aligned3.query_lst[l]
            current_ref3 = aligned3.ref_lst[l]
        elif l == len(aligned3.query_lst):
            current_seq3 = (-1,-1)


        if i == current_seq1[0]:
            is_mapped1 = 1
        if i == current_seq2[0]:
            is_mapped2 = 1
        if i == current_seq3[0]:
            is_mapped3 = 1
        #base mapped by all three aligners
        if is_mapped1 == 1 and is_mapped2 == 1 and is_mapped3 == 1:
            #all mapped to same position
            if ref1 == ref2 and ref2 == ref3 and current_ref1[0] == current_ref2[0] and current_ref2[0] == current_ref3[0]:
                same_pos_3 += 1
            #aligner1 and 2 map to same position
            elif ref1 == ref2 and current_ref1[0] == current_ref2[0]:
                same_pos_1_2_diff_3 += 1
            #aligner2 and 3 map to same position
            elif ref2 == ref3 and current_ref3[0] == current_ref2[0]:
                same_pos_2_3_diff_1 += 1
            #aligner1 and 3 map to same position
            elif ref1 == ref3 and current_ref1[0] == current_ref3[0]:
                same_pos_1_3_diff_2 += 1
            #no position are the same
            else:
                diff_pos_3 += 1
        #base mapped by only aligner 1 and 2
        elif is_mapped1 == 1 and is_mapped2 == 1 and is_mapped3 == 0:
            if ref1 == ref2 and current_ref1[0] == current_ref2[0]:
                same_pos_1_2_unmapped_3 += 1
            else:
                diff_pos_1_2_unmapped_3 += 1
        #base mapped by only aligner 1 and 3
        elif is_mapped1 == 1 and is_mapped3 == 1 and is_mapped2 == 0:
            if ref1 == ref3 and current_ref1[0] == current_ref3[0]:
                same_pos_1_3_unmapped_2 += 1
            else:
                diff_pos_1_3_unmapped_2 += 1
        #base mapped by only aligner 2 and 3
        elif is_mapped2 == 1 and is_mapped3 == 1 and is_mapped1 == 0:
            if ref2 == ref3 and current_ref2[0] == current_ref3[0]:
                same_pos_2_3_unmapped_1 += 1
            else:
                diff_pos_2_3_unmapped_1 += 1
        #mapped only by aligner 1
        elif is_mapped1 == 1:
            mapped_only_1 += 1
        #mapped only by aligner 2
        elif is_mapped2 == 1:
            mapped_only_2 += 1
        #mapped only by aligner 3
        elif is_mapped3 == 1:
            mapped_only_3 += 1
        #base not mapped by any aligner
        else:
            all_unmapped += 1

        #update positions in query and reference
        if is_mapped1 == 1:
            if i+1 == current_seq1[1]:
                j += 1
            else:
                aligned1.query_lst[j] = (i+1, current_seq1[1])
                aligned1.ref_lst[j] = (current_ref1[0]+1, current_ref1[1])

        #update positions in query and reference
        if is_mapped2 == 1:
            if i+1 == current_seq2[1]:
                k += 1
            else:
                aligned2.query_lst[k] = (i+1, current_seq2[1])
                aligned2.ref_lst[k] = (current_ref2[0]+1, current_ref2[1])
        #update positions in query and reference
        if is_mapped3 == 1:
            if i+1 == current_seq3[1]:
                l += 1
            else:
                aligned3.query_lst[l] = (i+1, current_seq3[1])
                aligned3.ref_lst[l] = (current_ref3[0]+1, current_ref3[1])

        i += 1



    return (same_pos_3, same_pos_1_2_diff_3, same_pos_2_3_diff_1,
        same_pos_1_3_diff_2, same_pos_1_2_unmapped_3, same_pos_2_3_unmapped_1,
        same_pos_1_3_unmapped_2, diff_pos_3, diff_pos_1_2_unmapped_3, diff_pos_2_3_unmapped_1,
        diff_pos_1_3_unmapped_2, mapped_only_1, mapped_only_2, mapped_only_3, all_unmapped)



def concordance_driver():
    assert len(sys.argv) == 2, "need to include mapping info file"
    mapping_info = sys.argv[1]

    """read the mapping_info csv file
    headers: row, qname, qlength, avearge quality,
    aligner1, ref1, ref_start1, CIGAR1,
    aligner2, ref2, ref_start2, CIGAR2,
    aligner3, ref3, ref_start3, CIGAR3
    """

    f1 = open(mapping_info, 'r')
    #output the concordance stats one query per line
    f2 = open("three_concordance.csv", 'w+')
    f2.write("qname, qlength, same_pos_3, same_pos_1_2_diff_3, same_pos_2_3_diff_1, \
    same_pos_1_3_diff_2, same_pos_1_2_unmapped_3, same_pos_2_3_unmapped_1,\
    same_pos_1_3_unmapped_2, diff_pos_3, diff_pos_1_2_unmapped_3, \
    diff_pos_2_3_unmapped_1, diff_pos_1_3_unmapped_2, mapped_only_1, mapped_only_2, mapped_only_3, all_unmapped\n")

    for line in f1.readlines()[1:]:

        query = line[:-1].split(',')
        qname = query[1][1:-1]
        qlength = int(query[2])
        aligned_by_1 = int(query[4])
        ref1 = query[5][1:-1]
        ref_start1 = query[6]
        CIGAR1 = query[7][1:-1]
        aligned_by_2 = int(query[8])
        ref2 = query[9][1:-1]
        ref_start2 = query[10]
        CIGAR2 = query[11][1:-1]
        aligned_by_3 = int(query[12])
        ref3 = query[13][1:-1]
        ref_start3 = query[14]
        CIGAR3 = query[15][1:-1]

        out_str = ""
        if aligned_by_1 == 1 and aligned_by_2 == 1 and aligned_by_3 == 1:
            alignment1 = process_CIGAR(CIGAR1, ref_start1)
            alignment2 = process_CIGAR(CIGAR2, ref_start2)
            alignment3 = process_CIGAR(CIGAR3, ref_start3)
            info = concordance_stat3(alignment1, alignment2, alignment3,qlength, ref1, ref2, ref3)
            info_str = str(info)[1:-1]
            out_str = format("%s,%d,%s\n"%(qname, qlength, info_str))
        elif aligned_by_1 == 1 and aligned_by_2 == 1:
            alignment1 = process_CIGAR(CIGAR1, ref_start1)
            alignment2 = process_CIGAR(CIGAR2, ref_start2)
            same_ref = 0
            if ref1 == ref2:
                same_ref = 1
            info = concordance_stat(alignment1, alignment2, qlength, same_ref)
            #unaligned_by_both,aligned_by_both, aligned_different_pos, aligned_1, aligned_2
            out_str = format("%s,%d,0,0,0,0,%d,0,0,0,%d,0,0,%d,%d,0,%d\n" \
                %(qname, qlength,info[1],info[2],info[3],info[4],info[0]))
        elif aligned_by_1 == 1 and aligned_by_3 == 1:
            alignment1 = process_CIGAR(CIGAR1, ref_start1)
            alignment3 = process_CIGAR(CIGAR3, ref_start3)
            same_ref = 0
            if ref1 == ref3:
                same_ref = 1
            info = concordance_stat(alignment1, alignment3, qlength, same_ref)
            out_str = format("%s,%d,0,0,0,0,0,0,%d,0,0,0,%d,%d,0,%d,%d\n" \
                %(qname, qlength, info[1], info[2],info[3], info[4], info[0]))
        elif aligned_by_2 == 1 and aligned_by_3 == 1:
            alignment2 = process_CIGAR(CIGAR2, ref_start2)
            alignment3 = process_CIGAR(CIGAR3, ref_start3)
            same_ref = 0
            if ref2 == ref3:
                same_ref = 1
            info = concordance_stat(alignment2, alignment3, qlength, same_ref)
            out_str = format("%s,%d,0,0,0,0,0,%d,0,0,0,%d,0,0,%d,%d,%d\n"
            %(qname, qlength,info[1], info[2],info[3], info[4], info[0]))
        elif aligned_by_1 == 1:
            aligned = aligned_base_count(CIGAR1)
            out_str = format("%s,%d,0,0,0,0,0,0,0,0,0,0,0,%d,0,0,%d\n"%(qname, qlength, aligned, qlength-aligned))
        elif aligned_by_2 == 1:
            aligned = aligned_base_count(CIGAR2)
            out_str = format("%s,%d,0,0,0,0,0,0,0,0,0,0,0,0,%d,0,%d\n"%(qname, qlength, aligned, qlength-aligned))
        elif aligned_by_3 == 1:
            aligned = aligned_base_count(CIGAR3)
            out_str = format("%s,%d,0,0,0,0,0,0,0,0,0,0,0,0,0,%d,%d\n"%(qname, qlength, aligned, qlength-aligned))
        else:
            out_str = format("%s,%d,0,0,0,0,0,0,0,0,0,0,0,0,0,0,%d\n"%(qname, qlength, qlength))
        f2.write(out_str)

    f1.close()
    f2.close()
    return

def main():
    #test the function
    qlength = 655
    pos1 = 12975116
    pos2 = 12975117
    pos3 = 128342167
    ref1 = "chrX"
    ref2 = "chrX"
    ref3 = "chr9"
    CIGAR1 = "30M1D3M1D5M1I2M1I6M1I4M1077N1M1D56M1I46M3D9M415N34M2D1M1I32M2D24M2D4M1D10M2D90M2D28M1D26M2D9M1D13M1D40M1D44M5D11M2I2M1D2M2D22M94S"
    CIGAR2 = "1S28M1D4M1D5M1I2M1I6M1I4M1077N1M1D55M1I48M2D1M1D7M415N34M2D1M1I34M1D1M1D25M3D11M2D46M1I3M1D23M1D4M1I12M2D29M1D24M2D13M1D13M1D38M1D48M5D5M2I2M1D5M1D3M1D16M5D3M2D5M1D4M1D3M1D5M3D1M73S"
    CIGAR3 = "4S25M1D2M1I6M1I1M2D1M1I6M1I5M1D55M1I46M1D1M2D41M1D2M1D1M1I29M1D3M1D22M3D16M2D47M1I1M1D25M1D2M1I11M1D1M1D8M3I2M1I15M1D13M252S"

    aligned1 = process_CIGAR(CIGAR1, pos1)
    aligned2 = process_CIGAR(CIGAR2, pos2)
    aligned3 = process_CIGAR(CIGAR3, pos3)
    #print(aligned1.query_lst)
    #print(aligned1.ref_lst)
    #print(aligned2.query_lst)
    #print(aligned2.ref_lst)
    #concordance_stat3(aligned1, aligned2, aligned3, qlength, ref1, ref2, ref3)
    concordance_driver()

if __name__ == "__main__":
    main()
