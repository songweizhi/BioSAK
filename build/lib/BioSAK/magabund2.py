import os
import glob
import argparse
from Bio import SeqIO
from itertools import groupby


magabund_usage = '''
=============================== magabund example commands ===============================

# example command 
BioSAK magabund -s all_bins.sam -m all_bins -x fa -o output.txt

# Note
The input sam file is obtained by mapping quality-filtered reads to the assemblies 
derived from these reads. 

=========================================================================================
'''


def Cdb_2_bin_cluster_file(Cdb_file, bin_cluster_file):
    cluster_to_bin_dict = {}
    obtained_clusters = set()
    for each_bin in open(Cdb_file):
        if not each_bin.startswith('genome,secondary_cluster'):
            each_bin_split = each_bin.strip().split(',')
            bin_id = each_bin_split[0]
            secondary_cluster = each_bin_split[1]
            obtained_clusters.add(secondary_cluster)
            if secondary_cluster not in cluster_to_bin_dict:
                cluster_to_bin_dict[secondary_cluster] = [bin_id]
            else:
                cluster_to_bin_dict[secondary_cluster].append(bin_id)

    obtained_clusters_list = sorted([i for i in obtained_clusters])

    bin_cluster_file_handle = open(bin_cluster_file, 'w')
    for j in obtained_clusters_list:
        bin_cluster_file_handle.write('cluster_%s\t%s\n' % (j, '\t'.join(cluster_to_bin_dict[j])))
    bin_cluster_file_handle.close()


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
        if each_part_cate in ['S', 's']:
            clipping_len += each_part_len

        # get mismatch_len
        if each_part_cate in {'I', 'i', 'X', 'x', 'D', 'd'}:
            mismatch_len += each_part_len

    aligned_pct  = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
    clipping_pct = float("{0:.2f}".format(clipping_len * 100 / (aligned_len + clipping_len)))
    mismatch_pct = float("{0:.2f}".format(mismatch_len * 100 / (aligned_len)))

    return aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct


# def get_ref_to_total_read_len_from_sam(input_sam_file, ctg2gnm_dict, read2len_dict):
#
#     ref_len_dict = dict()
#     reads_in_binned = set()
#     read_in_sam = set()
#     gnm_to_aligned_read_dict = dict()
#     for each_read in open(input_sam_file):
#         each_read_split = each_read.strip().split('\t')
#         if each_read.startswith('@'):
#             ref_id = ''
#             ref_len = 0
#             for each_element in each_read_split:
#                 if each_element.startswith('SN:'):
#                     ref_id = each_element[3:]
#                 if each_element.startswith('LN:'):
#                     ref_len = int(each_element[3:])
#             ref_len_dict[ref_id] = ref_len
#         else:
#             read_id   = each_read_split[0]
#             ref_id    = each_read_split[2]
#             ref_pos   = int(each_read_split[3])
#             cigar     = each_read_split[5]
#             cigar_splitted = cigar_splitter(cigar)
#             read_in_sam.add(read_id)
#             if ref_id == '*':
#                 gnm_id = 'unmapped'
#             else:
#                 gnm_id = ctg2gnm_dict.get(ref_id, 'unbinned')
#                 if gnm_id != 'unbinned':
#                     reads_in_binned.add(read_id)
#
#             to_ignore = False
#             both_ends_clipping = check_both_ends_clipping(cigar_splitted)
#             if both_ends_clipping is True:
#                 left_clip_len = int(cigar_splitted[0][:-1])
#                 left_clip_pct = left_clip_len*100/read2len_dict[read_id]
#                 right_clip_len = int(cigar_splitted[-1][:-1])
#                 right_clip_pct = right_clip_len*100/read2len_dict[read_id]
#                 if (left_clip_pct >=5) and (right_clip_pct >=5) and (left_clip_len >= 50) and (right_clip_len >=50):
#                     to_ignore = True
#
#             # check if clp in the middle
#             if ref_id == '*':
#                 clip_in_middle = False
#             else:
#                 aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
#                 ref_len = ref_len_dict[ref_id]
#                 clip_in_middle = False
#                 if ('S' in cigar) or ('s' in cigar):
#                     if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos > 200):
#                         clip_in_middle = True
#                     if (cigar_splitted[-1][-1] in ['S', 's']):
#                         if (ref_pos + aln_len - 1) - ref_len > 200:
#                             clip_in_middle = True
#
#             if (to_ignore is False) and (clip_in_middle is False):
#                 if gnm_id not in gnm_to_aligned_read_dict:
#                     gnm_to_aligned_read_dict[gnm_id] = {read_id}
#                 else:
#                     gnm_to_aligned_read_dict[gnm_id].add(read_id)
#
#     unbinned_updated = set()
#     for each_unbinned_read in gnm_to_aligned_read_dict.get('unbinned', set()):
#         if each_unbinned_read not in reads_in_binned:
#             unbinned_updated.add(each_unbinned_read)
#     gnm_to_aligned_read_dict['unbinned'] = unbinned_updated
#
#     ###############
#
#     gnm_to_aligned_read_len_dict = dict()
#     for each_gnm in gnm_to_aligned_read_dict:
#         gnm_to_aligned_read_len_dict[each_gnm] = 0
#         for each_read in gnm_to_aligned_read_dict[each_gnm]:
#             read_len = read2len_dict[each_read]
#             gnm_to_aligned_read_len_dict[each_gnm] += read_len
#
#     return gnm_to_aligned_read_len_dict, read_in_sam


def get_ref_to_total_read_len_from_sam(input_sam_file, ctg2gnm_dict, read2len_dict):

    read_in_sam = set()
    reads2gnm_dict = dict()
    for each_read in open(input_sam_file):
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
                gnm_id = ctg2gnm_dict.get(ref_id, 'unbinned')

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

    gnm_to_aligned_read_len_dict = dict()
    for each_gnm in gnm_to_aligned_read_dict:
        gnm_to_aligned_read_len_dict[each_gnm] = 0
        for each_read in gnm_to_aligned_read_dict[each_gnm]:
            read_len = read2len_dict[each_read]
            gnm_to_aligned_read_len_dict[each_gnm] += read_len

    return gnm_to_aligned_read_len_dict, read_in_sam


def magabund(args):

    ################################################# read in arguments ################################################

    sam_file      = args['s']
    bin_folder    = args['m']
    bin_ext       = args['x']
    output_file   = args['o']
    only_aligned  = args['a']
    # cluster_info  = args['g']
    # dRep_Cdb_file = args['Cdb']

    ############################################## define bin_cluster file #############################################

    # bin_cluster_file = ''
    # if (cluster_info is None) and (dRep_Cdb_file is None):
    #     bin_cluster_file = None
    #
    # elif (cluster_info is not None) and (dRep_Cdb_file is None):
    #     bin_cluster_file = cluster_info
    #
    # elif (cluster_info is None) and (dRep_Cdb_file is not None):
    #     Cdb_file_path, Cdb_file_basename, Cdb_file_extension = sep_path_basename_ext(dRep_Cdb_file)
    #     cluster_file_from_Cdb = '%s/%s_derived_cluster_file_%s%s' % (Cdb_file_path, Cdb_file_basename, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'), Cdb_file_extension)
    #     Cdb_2_bin_cluster_file(dRep_Cdb_file, cluster_file_from_Cdb)
    #     bin_cluster_file = cluster_file_from_Cdb
    # else:
    #     print(datetime.now().strftime(time_format) + 'cluster_info and dRep_Cdb are not compatible, please specify one only, program exited!')
    #     exit()

    ################################################ get ctg to bin dict ###############################################

    bin_file_re = '%s/*%s' % (bin_folder, bin_ext)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    if len(bin_file_list) == 0:
        print('No genome found, program exited!')
        exit()

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

    if len(ctgs_found_in_multiple_mags) > 0:
        print('The following %s sequences found in multiple genomes, program exited!' % len(ctgs_found_in_multiple_mags))
        print(','.join(ctgs_found_in_multiple_mags))
        exit()

    ############################################ get ctg to group list dict ############################################

    # print(datetime.now().strftime(time_format) + 'Get bin (cluster) to contig correlations')
    #
    # if bin_cluster_file is None:
    #     group_2_ctg_dict = bin_2_ctg_dict
    # else:
    #     # get group to bin dict
    #     group_2_bin_dict = {}
    #     for group in open(bin_cluster_file):
    #         group_split = group.strip().split('\t')
    #         group_2_bin_dict[group_split[0]] = group_split[1:]
    #
    #     group_2_ctg_dict = {}
    #     for bin_group in group_2_bin_dict:
    #         group_member_list = group_2_bin_dict[bin_group]
    #         group_2_ctg_dict[bin_group] = set()
    #         for genome_bin in group_member_list:
    #             genome_bin_ctg_list = bin_2_ctg_dict[genome_bin]
    #             for ctg in genome_bin_ctg_list:
    #                 group_2_ctg_dict[bin_group].add(ctg)

    ########################################### get_ref_to_read_num_from_sam ###########################################

    reads_file = 'assembly_graph_NoEu.gfa_reads_min0bp.fastq'

    read_len_dict = dict()
    for each_read in SeqIO.parse(reads_file, 'fastq'):
        read_len_dict[each_read.id] = len(each_read.seq)

    print('Get the total length of reads mapped to each reference sequence in sam file')

    gnm_to_aligned_read_total_len_dict, reads_in_sam = get_ref_to_total_read_len_from_sam(sam_file, ctg_2_bin_dict, read_len_dict)

    reads_not_in_sam_total_len = 0
    for each_read in read_len_dict:
        read_len = read_len_dict[each_read]
        if each_read not in reads_in_sam:
            reads_not_in_sam_total_len += read_len
    print('reads_not_in_sam_total_len\t%s' % reads_not_in_sam_total_len)

    if 'unmapped' not in gnm_to_aligned_read_total_len_dict:
        gnm_to_aligned_read_total_len_dict['unmapped'] = reads_not_in_sam_total_len
    else:
        gnm_to_aligned_read_total_len_dict['unmapped'] += reads_not_in_sam_total_len


    total_read_len    = sum(gnm_to_aligned_read_total_len_dict.values())
    total_gnm_size_bp = sum(bin_size_dict.values())
    gnm_list          = sorted(list(gnm_to_aligned_read_total_len_dict.keys()))

    # get total_mag_depth
    total_mag_depth = 0
    for each_gnm in gnm_list:
        if each_gnm not in ['unbinned', 'unmapped']:
            gnm_size_bp     = bin_size_dict[each_gnm]
            gnm_reads_len   = gnm_to_aligned_read_total_len_dict[each_gnm]
            gnm_reads_depth = float("{0:.2f}".format(gnm_reads_len/gnm_size_bp))
            total_mag_depth += gnm_reads_depth
    total_mag_depth = float("{0:.2f}".format(total_mag_depth))

    # write out table
    gnm_reads_len_pct_total = 0
    gnm_reads_depth_norm_total = 0
    output_file_handle = open(output_file, 'w')
    output_file_handle.write('Genome\tGenome_size(bp)\tMapped_reads_length(bp)\tMapped_reads_percent\tDepth(x)\tDepth/%s\n' % total_mag_depth)
    for each_gnm in gnm_list:
        if each_gnm not in ['unbinned', 'unmapped']:
            gnm_size_bp                 = bin_size_dict[each_gnm]
            gnm_reads_len               = gnm_to_aligned_read_total_len_dict[each_gnm]
            gnm_reads_len_pct           = float("{0:.2f}".format(gnm_reads_len*100/total_read_len))
            gnm_reads_depth             = float("{0:.2f}".format(gnm_reads_len/gnm_size_bp))
            gnm_reads_depth_norm        = float("{0:.2f}".format(gnm_reads_depth*100/total_mag_depth))
            gnm_reads_len_pct_total     += gnm_reads_len_pct
            gnm_reads_depth_norm_total  += gnm_reads_depth_norm
            output_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_gnm, gnm_size_bp, gnm_reads_len, gnm_reads_len_pct, gnm_reads_depth, gnm_reads_depth_norm))

    if 'unbinned' in gnm_to_aligned_read_total_len_dict:
        gnm_reads_len           = gnm_to_aligned_read_total_len_dict['unbinned']
        gnm_reads_len_pct       = float("{0:.2f}".format(gnm_reads_len*100/total_read_len))
        gnm_reads_len_pct_total += gnm_reads_len_pct
        output_file_handle.write('%s\tNA\t%s\t%s\tNA\tNA\n' % ('unbinned', gnm_reads_len, gnm_reads_len_pct))

    if 'unmapped' in gnm_to_aligned_read_total_len_dict:
        gnm_reads_len           = gnm_to_aligned_read_total_len_dict['unmapped']
        gnm_reads_len_pct       = float("{0:.2f}".format(gnm_reads_len*100/total_read_len))
        gnm_reads_len_pct_total += gnm_reads_len_pct
        output_file_handle.write('%s\tNA\t%s\t%s\tNA\tNA\n' % ('unmapped', gnm_reads_len, gnm_reads_len_pct))

    output_file_handle.write('Total\t%s\t%s\t%s\t%s\t%s\n' % (total_gnm_size_bp, total_read_len, gnm_reads_len_pct_total, float("{0:.2f}".format(total_mag_depth)), gnm_reads_depth_norm_total))
    output_file_handle.close()

    # final report
    print('Done!')


if __name__ == "__main__":

    magabund_parser = argparse.ArgumentParser()
    magabund_parser.add_argument('-s',     required=True,                        help='input sam file')
    magabund_parser.add_argument('-m',     required=True,                        help='MAG folder')
    magabund_parser.add_argument('-x',     required=True, default='fasta',       help='MAG file extension, default: fasta')
    magabund_parser.add_argument('-o',     required=True,                        help='output table')
    magabund_parser.add_argument('-a',     required=False, action='store_true',  help='based on aligned length, rather than read length')
    #magabund_parser.add_argument('-g',       required=False, default=None,         help='bin grouping info')
    #magabund_parser.add_argument('-Cdb',     required=False, default=None,         help='cluster info from dRep (Cdb.csv)')
    args = vars(magabund_parser.parse_args())
    magabund(args)

'''
# get abundance of customized bin clusters
BioSAK get_bin_abundance -sam all_bins.sam -bin all_bins -x fa -o output.txt -g bin_clusters.txt

# get abundance of dRep produced bin clusters
BioSAK get_bin_abundance -sam all_bins.sam -bin all_bins -x fa -o output.txt -Cdb Cdb.csv

cd /Users/songweizhi/Desktop/get_MAG_abundance_wd
python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/magabund.py -sam three_drep_mags.sam -bin mags -x fa -o test.txt

# How it works
The get_bin_abundance module first gets the number of reads mapped to each reference sequence in the provided sam file. 
Then calculate the number/percentage of reads mapped to the sequences in each bin (cluster).

# Format of customized grouping file (tab-separated, with the first column as group id, followed by a list of bins from it)
cluster_1	Ecklonia_bin_1.fa	Delisea_bin_1.fa
cluster_2	Ecklonia_bin_2.fa	Delisea_bin_5.fa	Amphiroa_bin_1.fa
cluster_3	Amphiroa_bin_3.fa

'''
