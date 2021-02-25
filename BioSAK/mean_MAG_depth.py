import os
import glob
import argparse
from Bio import SeqIO
from datetime import datetime


get_bin_abundance_usage = '''
========================== get_bin_abundance example commands ==========================

# Example command
BioSAK get_bin_abundance -d ctg_lt2500_depth.txt -b bin_files -x fasta -p Refined_bins

# How it works

========================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_bin_abundance(args):

    # read in arguments
    metabat_depth_file =  args['d']
    bin_folder =          args['b']
    bin_file_extension =  args['x']
    output_file =         args['p']

    # check input file
    if bin_folder[-1] == '/':
      bin_folder = bin_folder[:-1]

    # get bin file list
    bin_file_re= '%s/*.%s' % (bin_folder, bin_file_extension)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    # check whether input bin files detected
    if len(bin_file_list) == 0:
        print('%s %s' % ((datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')), 'No bin file detected, program exit!'))
        exit()

    # read in depth info
    # get ctg_len_dict and ctg_depth_dict
    ctg_len_dict = {}
    ctg_depth_dict = {}
    for each_line in open(metabat_depth_file):
        if not each_line.startswith('contigName	contigLen	totalAvgDepth'):
            each_line_split = each_line.strip().split('\t')
            ctg_id = each_line_split[0]
            ctg_len = int(each_line_split[1])
            ctg_depth = float(each_line_split[2])
            ctg_len_dict[ctg_id] = ctg_len
            ctg_depth_dict[ctg_id] = ctg_depth

    file_re = '%s/*.%s' % (bin_folder, bin_file_extension)
    file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]

    gnm_mean_depth_dict = {}
    gnm_mean_depth_file_handle = open(output_file, 'w')
    for genome in file_list:
        pwd_genome = '%s/%s' % (bin_folder, genome)
        gnm_total_len = 0
        gnm_total_depth = 0
        for each_ctg in SeqIO.parse(pwd_genome, 'fasta'):
            gnm_total_len += len(each_ctg.seq)
            gnm_total_depth += (len(each_ctg.seq) * ctg_depth_dict[each_ctg.id])
        gnm_mean_depth = gnm_total_depth / gnm_total_len
        gnm_mean_depth = float("{0:.2f}".format(gnm_mean_depth))
        gnm_mean_depth_file_handle.write('%s\t%s\t%s\n' % (genome, gnm_total_len, gnm_mean_depth))
        gnm_mean_depth_dict[genome] = gnm_mean_depth
    gnm_mean_depth_file_handle.close()

    # report
    print('%s %s' % ((datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')), 'Done!'))


if __name__ == "__main__":

    get_bin_abundance_parser = argparse.ArgumentParser()

    # Annotation modules
    get_bin_abundance_parser.add_argument('-d', required=True,                    help='MetaBAT produced depth file')
    get_bin_abundance_parser.add_argument('-b', required=True,                    help='bin folder')
    get_bin_abundance_parser.add_argument('-x', required=False, default='fasta',  help='file extension')
    get_bin_abundance_parser.add_argument('-p', required=True,                    help='output file')

    args = vars(get_bin_abundance_parser.parse_args())

    get_bin_abundance(args)

