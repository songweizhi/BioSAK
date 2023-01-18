import os
import glob
import argparse
from Bio import SeqIO
from itertools import groupby

def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if (each_element.isalpha() is True) or (each_element == '='):
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def get_cigar_stats(cigar_splitted):

    # aligned_len: M I X =
    # clipping_len: S
    # mismatch_len: X I D
    # mismatch_pct = mismatch_len / aligned_len
    # aligned_pct  = aligned_len  / (aligned_len + clipping_len)
    # clipping_pct = clipping_len / (aligned_len + clipping_len)

    aligned_len = 0
    clipping_len = 0
    mismatch_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len

        # get clipping_len
        if each_part_cate in ['S', 's', 'H', 'h']:
            clipping_len += each_part_len

        # get mismatch_len
        if each_part_cate in {'I', 'i', 'X', 'x', 'D', 'd'}:
            mismatch_len += each_part_len

    aligned_pct  = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
    clipping_pct = float("{0:.2f}".format(clipping_len * 100 / (aligned_len + clipping_len)))
    mismatch_pct = float("{0:.2f}".format(mismatch_len * 100 / (aligned_len)))

    return aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct



sam_file    = '/Users/songweizhi/Desktop/test/map_to.fa.no.secondary.sam'
sam_file    = '/Users/songweizhi/Desktop/test/test.sam'
bin_folder  = '/Users/songweizhi/Desktop/test/metawrap_30_90_bins_with_completed_renamed'
bin_ext     = 'fa'
output_file = '/Users/songweizhi/Desktop/test/abund.txt'

bin_file_re = '%s/*%s' % (bin_folder, bin_ext)
bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

bin_size_dict = dict()
ctg_2_bin_dict = dict()
counted_ctg_set = set()
ctgs_found_in_multiple_mags = set()
for each_bin in bin_file_list:
    pwd_each_bin = '%s/%s' % (bin_folder, each_bin)
    bin_size_dict[each_bin] = 0
    for seq in SeqIO.parse(pwd_each_bin, 'fasta'):
        ctg_2_bin_dict[seq.id] = each_bin
        bin_size_dict[each_bin] += len(seq.seq)
        if seq.id not in counted_ctg_set:
            counted_ctg_set.add(seq.id)
        else:
            ctgs_found_in_multiple_mags.add(seq.id)


read_in_sam = set()
reads2gnm_dict = dict()
for each_read in open(sam_file):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        ref_id = each_read_split[2]
        ref_pos = int(each_read_split[3])
        cigar = each_read_split[5]
        read_in_sam.add(read_id)
        if ref_id == '*':
            gnm_id = 'unmapped'
        else:
            gnm_id = ctg_2_bin_dict.get(ref_id, 'unbinned')

        if read_id not in reads2gnm_dict:
            reads2gnm_dict[read_id] = dict()
            reads2gnm_dict[read_id][gnm_id] = [cigar]
        else:
            if gnm_id not in reads2gnm_dict[read_id]:
                reads2gnm_dict[read_id][gnm_id] = [cigar]
            else:
                reads2gnm_dict[read_id][gnm_id].append(cigar)

gnm_to_aligned_read_dict = dict()
for each_read in reads2gnm_dict:
    mapped_gnms = reads2gnm_dict[each_read]

    if len(mapped_gnms) == 1:
        for mapped_gnm in mapped_gnms:
            if mapped_gnm not in gnm_to_aligned_read_dict:
                gnm_to_aligned_read_dict[mapped_gnm] = {each_read}
            else:
                gnm_to_aligned_read_dict[mapped_gnm].add(each_read)
    elif (len(mapped_gnms) == 2) and ('unbinned' in mapped_gnms):
        for mapped_gnm in mapped_gnms:
            if mapped_gnm != 'unbinned':
                if mapped_gnm not in gnm_to_aligned_read_dict:
                    gnm_to_aligned_read_dict[mapped_gnm] = {each_read}
                else:
                    gnm_to_aligned_read_dict[mapped_gnm].add(each_read)
    else:
        gnm_with_longest_aln_len = ''
        longest_aln_len = 0
        for mapped_gnm in mapped_gnms:
            if mapped_gnm != 'unbinned':
                for each_cigar in mapped_gnms[mapped_gnm]:
                    cigar_splitted = cigar_splitter(each_cigar)
                    aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                    if aln_len > longest_aln_len:
                        longest_aln_len = aln_len
                        gnm_with_longest_aln_len = mapped_gnm

        if gnm_with_longest_aln_len not in gnm_to_aligned_read_dict:
            gnm_to_aligned_read_dict[gnm_with_longest_aln_len] = {each_read}
        else:
            gnm_to_aligned_read_dict[gnm_with_longest_aln_len].add(each_read)










