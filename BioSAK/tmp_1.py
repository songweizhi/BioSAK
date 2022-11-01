from itertools import groupby


# sam_file = '/Users/songweizhi/Desktop/assembly.fasta.sam'
# read_id_txt = '/Users/songweizhi/Desktop/reads.txt'
# ref_id_txt = '/Users/songweizhi/Desktop/refs.txt'
#
# sam_file = 'assembly.fasta.sam'
# read_id_txt = 'reads.txt'
# ref_id_txt = 'refs.txt'
#
# interested_read_set = set()
# for each_read in open(read_id_txt):
#     interested_read_set.add(each_read.strip())
#
#
# ref_set = set()
# for each_aln in open(sam_file):
#     each_aln_split = each_aln.strip().split('\t')
#     read_id = each_aln_split[0]
#     ref_id = each_aln_split[2]
#     if read_id in interested_read_set:
#         ref_set.add(ref_id)
#
#
# ref_id_txt_handle = open(ref_id_txt, 'w')
# for each_ref in ref_set:
#     ref_id_txt_handle.write(each_ref + '\n')
# ref_id_txt_handle.close()


def cigar_to_read_len(cigar_string):
    # Given a CIGAR string, return the number of bases consumed from the query sequence
    read_consuming_ops = ("M", "I", "S", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            result += length
    return result



print(cigar_to_read_len('23=7I2S'))