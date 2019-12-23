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
    output_prefix =       args['p']

    # check input file
    if bin_folder[-1] == '/':
      bin_folder = bin_folder[:-1]

    # define output abundance file name
    output_abundance_txt = '%s_relative_abundance.txt' % output_prefix

    # get bin file list
    bin_file_re= '%s/*.%s' % (bin_folder, bin_file_extension)
    bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_file_re)]

    # check whether input bin files detected
    if len(bin_file_list) == 0:
        print('%s %s' % ((datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')), 'No bin file detected, program exit!'))
        exit()

    # read in depth info
    ctg_to_length_dict = {}
    ctg_to_depth_dict = {}
    ctg_to_weighted_depth_dict = {}
    line_num = 0
    ctg_total_depth_weighted = 0
    total_ctg_list = []
    for ctg_depth in open(metabat_depth_file):

        if line_num > 0:
            ctg_depth_split = ctg_depth.strip().split()
            ctg_id = ctg_depth_split[0]
            ctg_len = int(ctg_depth_split[1])
            ctg_depth = float(ctg_depth_split[2])
            ctg_depth_weighted = ctg_len*ctg_depth

            ctg_to_depth_dict[ctg_id] = ctg_depth
            ctg_to_length_dict[ctg_id] = ctg_len
            ctg_to_weighted_depth_dict[ctg_id] = ctg_depth_weighted
            ctg_total_depth_weighted += ctg_depth_weighted
            total_ctg_list.append(ctg_id)

        line_num += 1

    # write out depth info
    output_txt_handle = open(output_abundance_txt, 'w')
    output_txt_handle.write('Bin_id\tRelative_abundance\n')

    binned_ctg_list = set()
    for genome_bin in sorted(bin_file_list):

        pwd_genome_bin = '%s/%s' % (bin_folder, genome_bin)
        genome_bin_path, genome_bin_basename, genome_bin_extension = sep_path_basename_ext(pwd_genome_bin)

        current_bin_weighted_depth = 0
        for ctg in SeqIO.parse(pwd_genome_bin, 'fasta'):
            ctg_id = str(ctg.id)
            ctg_weighted_depth = ctg_to_weighted_depth_dict[ctg_id]
            current_bin_weighted_depth += ctg_weighted_depth
            binned_ctg_list.add(ctg_id)

        current_bin_weighted_depth_norm = float("{0:.2f}".format(current_bin_weighted_depth*100/ctg_total_depth_weighted))

        output_txt_handle.write('%s\t%s\n' % (genome_bin_basename, current_bin_weighted_depth_norm))

    # get unbinned contig list
    unbinned_ctg_list = []
    for ctg in total_ctg_list:
        if ctg not in binned_ctg_list:
            unbinned_ctg_list.append(ctg)

    # get unbinned contig total weighted depth
    unbinned_ctg_total_weighted_depth = 0
    for unbinned_ctg in unbinned_ctg_list:
        unbinned_ctg_weighted_depth = ctg_to_weighted_depth_dict[unbinned_ctg]
        unbinned_ctg_total_weighted_depth += unbinned_ctg_weighted_depth

    unbinned_ctg_total_weighted_depth_norm = unbinned_ctg_total_weighted_depth*100/ctg_total_depth_weighted
    unbinned_ctg_total_weighted_depth_norm = float("{0:.2f}".format(unbinned_ctg_total_weighted_depth_norm))

    output_txt_handle.write('%s\t%s\n' % ('unbinned', unbinned_ctg_total_weighted_depth_norm))
    output_txt_handle.close()

    # report
    print('%s %s' % ((datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')), 'Bin abundance exported to %s' % output_abundance_txt))
    print('%s %s' % ((datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')), 'Done!'))


if __name__ == "__main__":

    get_bin_abundance_parser = argparse.ArgumentParser()

    # Annotation modules
    get_bin_abundance_parser.add_argument('-d', required=True,                    help='MetaBAT produced depth file')
    get_bin_abundance_parser.add_argument('-b', required=True,                    help='bin folder')
    get_bin_abundance_parser.add_argument('-x', required=False, default='fasta',  help='file extension')
    get_bin_abundance_parser.add_argument('-p', required=False, default='OUTPUT', help='output prefix')

    args = vars(get_bin_abundance_parser.parse_args())

    get_bin_abundance(args)

