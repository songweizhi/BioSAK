import os
import glob
import argparse
from Bio import SeqIO


mean_MAG_cov_usage = '''
==================== mean_MAG_cov example commands ====================

# Example command
BioSAK mean_MAG_cov -d ctg_lt2500_depth.txt -b bin_folder -x fasta

Format of contig depth file: tab separated, no header!!!
ctg_1   8
ctg_2   26.1

=======================================================================
'''


def mean_MAG_cov(args):

    # read in arguments
    metabat_depth_file = args['d']
    bin_folder         = args['b']
    bin_file_extension = args['x']

    # check input file
    if bin_folder[-1] == '/':
      bin_folder = bin_folder[:-1]

    # get bin file list
    bin_file_re = '%s/*.%s' % (bin_folder, bin_file_extension)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    # check whether input bin files detected
    if len(bin_file_list) == 0:
        print('No input file detected, program exit!')
        exit()

    # read in depth info
    ctg_depth_dict = {}
    for each_line in open(metabat_depth_file):
        each_line_split = each_line.strip().split('\t')
        ctg_id = each_line_split[0]
        ctg_depth = float(each_line_split[1])
        ctg_depth_dict[ctg_id] = ctg_depth

    file_re = '%s/*.%s' % (bin_folder, bin_file_extension)
    file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]

    gnm_mean_depth_dict = {}
    print('MAG\tLength(bp)\tDepth')
    for genome in file_list:
        pwd_genome = '%s/%s' % (bin_folder, genome)
        gnm_total_len = 0
        gnm_total_depth = 0
        for each_ctg in SeqIO.parse(pwd_genome, 'fasta'):
            gnm_total_len += len(each_ctg.seq)
            gnm_total_depth += (len(each_ctg.seq) * ctg_depth_dict[each_ctg.id])
        gnm_mean_depth = gnm_total_depth / gnm_total_len
        gnm_mean_depth = float("{0:.2f}".format(gnm_mean_depth))
        print('%s\t%s\t%s' % (genome, gnm_total_len, gnm_mean_depth))
        gnm_mean_depth_dict[genome] = gnm_mean_depth


if __name__ == "__main__":
    mean_MAG_cov_parser = argparse.ArgumentParser()
    mean_MAG_cov_parser.add_argument('-d', required=True, help='MetaBAT produced depth file')
    mean_MAG_cov_parser.add_argument('-b', required=True, help='bin folder')
    mean_MAG_cov_parser.add_argument('-x', required=True, help='file extension')
    args = vars(mean_MAG_cov_parser.parse_args())
    mean_MAG_cov(args)

