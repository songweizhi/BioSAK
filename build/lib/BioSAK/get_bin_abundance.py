import os
import glob
import argparse
from Bio import SeqIO
from datetime import datetime
from BioSAK.global_functions import time_format
from BioSAK.global_functions import sep_path_basename_ext


get_bin_abundance_usage = '''
================================ get_bin_abundance example commands ================================

# get abundance of individual bins 
BioSAK get_bin_abundance -sam all_bins.sam -bin all_bins -x fa -o abundance.txt

# get abundance of customized bin clusters
BioSAK get_bin_abundance -sam all_bins.sam -bin all_bins -x fa -o abundance.txt -g bin_clusters.txt

# get abundance of dRep produced bin clusters
BioSAK get_bin_abundance -sam all_bins.sam -bin all_bins -x fa -o abundance.txt -Cdb Cdb.csv

# How it works
The get_bin_abundance module first gets the number of reads mapped to each reference sequence in the provided sam file. 
Then calculate the number/percentage of reads mapped to the sequences in each bin (cluster).

# Format of customized grouping file (tab-separated, with the first column as group id, followed by a list of bins from it)
cluster_1	Ecklonia_bin_1.fa	Delisea_bin_1.fa
cluster_2	Ecklonia_bin_2.fa	Delisea_bin_5.fa	Amphiroa_bin_1.fa
cluster_3	Amphiroa_bin_3.fa

# Note!!!
The input sam file is obtained by mapping sequencing reads from a sample to the combined file of all bins in ‘-bin’ folder. 
Please make sure contig ids are UNIQUE across all bins in ‘-bin’ folder.

====================================================================================================
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


def get_ref_to_read_num_from_sam(input_sam_file, output_stats_file):

    output_stats_file_path, output_stats_file_basename, output_stats_file_extension = sep_path_basename_ext(output_stats_file)
    output_stats_tmp = '%s/%s_tmp%s' % (output_stats_file_path, output_stats_file_basename, output_stats_file_extension)

    # Store reads num in dict
    ref2read_num_dict = {}
    for each_read in open(input_sam_file):
        if not each_read.startswith('@'):
            ref_id = each_read.strip().split('\t')[2]

            if ref_id not in ref2read_num_dict:
                ref2read_num_dict[ref_id] = 1
            else:
                ref2read_num_dict[ref_id] += 1

    # Write reads num to file
    stat_file_unsorted_handle = open(output_stats_tmp, 'w')
    for each_ref in ref2read_num_dict:
        stat_file_unsorted_handle.write('%s\t%s\n' % (each_ref, ref2read_num_dict[each_ref]))
    stat_file_unsorted_handle.close()

    # sort output file
    os.system('cat %s | sort > %s' % (output_stats_tmp, output_stats_file))

    # remove tmp file
    os.system('rm %s' % output_stats_tmp)


def get_bin_abundance(args):

    ################################################# read in arguments ################################################

    sam_file      = args['sam']
    bin_folder    = args['bin']
    bin_ext       = args['x']
    output_file   = args['o']
    cluster_info  = args['g']
    dRep_Cdb_file = args['Cdb']

    ############################################## define bin_cluster file #############################################

    bin_cluster_file = ''
    if (cluster_info is None) and (dRep_Cdb_file is None):
        bin_cluster_file = None

    elif (cluster_info is not None) and (dRep_Cdb_file is None):
        bin_cluster_file = cluster_info

    elif (cluster_info is None) and (dRep_Cdb_file is not None):

        Cdb_file_path, Cdb_file_basename, Cdb_file_extension = sep_path_basename_ext(dRep_Cdb_file)
        cluster_file_from_Cdb = '%s/%s_derived_cluster_file_%s%s' % (Cdb_file_path, Cdb_file_basename, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'), Cdb_file_extension)
        Cdb_2_bin_cluster_file(dRep_Cdb_file, cluster_file_from_Cdb)
        bin_cluster_file = cluster_file_from_Cdb

    else:
        print(datetime.now().strftime(time_format) + 'cluster_info and dRep_Cdb are not compatible, please specify one only, program exited!')
        exit()

    ################################################ get bin to ctg dict ###############################################

    bin_file_re = '%s/*%s' % (bin_folder, bin_ext)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    if len(bin_file_list) == 0:
        print(datetime.now().strftime(time_format) + 'No bin file found, program exited!')
        exit()

    bin_2_ctg_dict = {}
    for each_bin in bin_file_list:
        pwd_each_bin = '%s/%s' % (bin_folder, each_bin)
        bin_2_ctg_dict[each_bin] = set()
        for seq in SeqIO.parse(pwd_each_bin, 'fasta'):
            bin_2_ctg_dict[each_bin].add(seq.id)

    ############################################ get group to ctg list dict ############################################

    print(datetime.now().strftime(time_format) + 'Get bin (cluster) to contig correlations')

    if bin_cluster_file is None:
        group_2_ctg_dict = bin_2_ctg_dict
    else:
        # get group to bin dict
        group_2_bin_dict = {}
        for group in open(bin_cluster_file):
            group_split = group.strip().split('\t')
            group_2_bin_dict[group_split[0]] = group_split[1:]

        group_2_ctg_dict = {}
        for bin_group in group_2_bin_dict:
            group_member_list = group_2_bin_dict[bin_group]
            group_2_ctg_dict[bin_group] = set()
            for genome_bin in group_member_list:
                genome_bin_ctg_list = bin_2_ctg_dict[genome_bin]
                for ctg in genome_bin_ctg_list:
                    group_2_ctg_dict[bin_group].add(ctg)

    ########################################### get_ref_to_read_num_from_sam ###########################################

    print(datetime.now().strftime(time_format) + 'Get the number of reads mapped to each reference sequence in sam file')

    sam_file_path, sam_file_basename, sam_file_extension = sep_path_basename_ext(sam_file)
    ref_to_read_num_file = '%s/%s_ref_to_read_num_%s.txt' % (sam_file_path, sam_file_basename, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))

    get_ref_to_read_num_from_sam(sam_file, ref_to_read_num_file)

    ########################################### read in ref_to_read_num_file ###########################################

    mapped_reads_num = 0
    ref_to_read_num_dict = {}
    for each_ctg in open(ref_to_read_num_file):
        each_ctg_split = each_ctg.strip().split('\t')
        ctg_id = each_ctg_split[0]
        read_num = int(each_ctg_split[1])
        ref_to_read_num_dict[ctg_id] = read_num
        if ctg_id != '*':
            mapped_reads_num += read_num

    ###################################### get the number of reads in each group #######################################

    print(datetime.now().strftime(time_format) + 'Get the number of reads mapped to each bin (cluster)')

    group_to_read_num_dict = {}
    for group in group_2_ctg_dict:
        group_ctg_list = group_2_ctg_dict[group]
        group_to_read_num_dict[group] = 0
        for ctg in group_ctg_list:
            group_to_read_num_dict[group] += ref_to_read_num_dict.get(ctg, 0)

    output_file_handle = open(output_file, 'w')
    output_file_handle.write('cluster\tread_num\tread_pct\n')
    for group in group_to_read_num_dict:
        group_read_num = group_to_read_num_dict[group]
        group_read_pct = float("{0:.2f}".format(group_read_num*100/mapped_reads_num))
        output_file_handle.write('%s\t%s\t%s\n' % (group, group_read_num, group_read_pct))
    output_file_handle.close()

    ################################################## final report ####################################################

    # delete tmp files
    os.system('rm %s' % ref_to_read_num_file)
    if (cluster_info is None) and (dRep_Cdb_file is not None):
        os.system('rm %s' % bin_cluster_file)

    # final report
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == "__main__":

    get_bin_abundance_parser = argparse.ArgumentParser()

    # Annotation modules
    get_bin_abundance_parser.add_argument('-sam',     required=True,                  help='input sam file')
    get_bin_abundance_parser.add_argument('-bin',     required=True,                  help='bin folder')
    get_bin_abundance_parser.add_argument('-x',       required=True, default='fasta', help='bin file extension, default: fasta')
    get_bin_abundance_parser.add_argument('-o',       required=True,                  help='output abundance file')
    get_bin_abundance_parser.add_argument('-g',       required=False, default=None,   help='bin grouping info')
    get_bin_abundance_parser.add_argument('-Cdb',     required=False, default=None,   help='cluster info from dRep (Cdb.csv)')

    args = vars(get_bin_abundance_parser.parse_args())

    get_bin_abundance(args)
