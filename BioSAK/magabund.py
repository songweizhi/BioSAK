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

def cigar_to_aln_len(cigar_string):

    # Given a CIGAR string, return the number of bases consumed from the query sequence
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in ("M", "=", "X"):
            result += length
    return result


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


def get_ref_to_total_read_len_from_sam(input_sam_file, by_aligned):

    ctg_to_aligned_read_len_dict = {}
    for each_read in open(input_sam_file):
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            ref_id = each_read_split[2]
            cigar_str = each_read_split[5]
            read_seq = each_read_split[9]
            read_len = len(read_seq)
            if read_seq == '*':

                if by_aligned is True:
                    read_len = cigar_to_aln_len(cigar_str)
                else:
                    read_len = cigar_to_read_len(cigar_str)

            if ref_id not in ctg_to_aligned_read_len_dict:
                ctg_to_aligned_read_len_dict[ref_id] = read_len
            else:
                ctg_to_aligned_read_len_dict[ref_id] += read_len

    return ctg_to_aligned_read_len_dict


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

    print('Get the total length of reads mapped to each reference sequence in sam file')

    ctg_to_aligned_read_total_len_dict = get_ref_to_total_read_len_from_sam(sam_file, only_aligned)

    gnm_to_aligned_read_total_len_dict = dict()
    for each_ctg in ctg_to_aligned_read_total_len_dict:
        aligned_read_total_len = ctg_to_aligned_read_total_len_dict[each_ctg]
        if each_ctg == '*':
            gnm_to_aligned_read_total_len_dict['unmapped'] = aligned_read_total_len
        else:
            ctg_gnm = ctg_2_bin_dict.get(each_ctg, 'unbinned')
            if ctg_gnm not in gnm_to_aligned_read_total_len_dict:
                gnm_to_aligned_read_total_len_dict[ctg_gnm] = aligned_read_total_len
            else:
                gnm_to_aligned_read_total_len_dict[ctg_gnm] += aligned_read_total_len

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
    output_file_handle.write('Genome\tGenome_size(bp)\tTotal_read_length(bp)\tTotal_read_length(pct)\tDepth(x)\tDepth/%s(Total_MAG_Depth)\n' % total_mag_depth)
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
