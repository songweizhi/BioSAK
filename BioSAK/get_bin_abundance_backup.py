import os
import glob
import argparse
from Bio import SeqIO
from datetime import datetime


get_bin_abundance_usage = '''
==================== get_bin_abundance example commands ====================

# Example command
BioSAK get_bin_abundance -d ctg_lt2500_depth.txt -b bin_files -x fasta -p Refined_bins

# Software dependencies:
module load python/3.7.3

==========================================================================
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

    # store info in dicts
    bin_list = []
    bin_to_ctg_dict = {}
    ctg_to_bin_dict = {}
    ctg_to_len_dict = {}
    bin_to_size_dict = {}
    for genome_bin in bin_file_list:

        pwd_genome_bin = '%s/%s' % (bin_folder, genome_bin)
        genome_bin_path, genome_bin_basename, genome_bin_extension = sep_path_basename_ext(pwd_genome_bin)

        current_bin_ctg_list = []
        current_bin_size = 0
        for contig in SeqIO.parse(pwd_genome_bin, 'fasta'):

            contig_id = str(contig.id)
            contig_len = len(str(contig.seq))

            current_bin_size += contig_len
            ctg_to_bin_dict[contig_id] = genome_bin_basename
            ctg_to_len_dict[contig_id] = contig_len
            current_bin_ctg_list.append(contig_id)

        bin_list.append(genome_bin_basename)
        bin_to_ctg_dict[genome_bin_basename] = sorted(current_bin_ctg_list)
        bin_to_size_dict[genome_bin_basename] = current_bin_size

    # read in depth info
    sample_name_list = []
    sample_num = 0
    ctg_to_depth_dict = {}
    line_num = 0
    for ctg_depth in open(metabat_depth_file):

        ctg_depth_split = ctg_depth.strip().split()

        # get sample_name_list
        if line_num == 0:
            sample_depth_header = ctg_depth_split[3:]
            sample_index = 0
            while sample_index < len(sample_depth_header):
                sample_name_list.append(sample_depth_header[sample_index])
                sample_index += 2
        else:
            ctg_id = ctg_depth_split[0]
            depth_list = ctg_depth_split[3:]
            sample_num = int(len(depth_list)/2)

            current_ctg_depth_list = []
            n = 0
            while n < sample_num * 2:
                current_sample_depth = float(depth_list[n])
                current_ctg_depth_list.append(current_sample_depth)
                n += 2

            #print('%s\t%s' % (ctg_id, current_ctg_depth_list))

            ctg_to_depth_dict[ctg_id] = current_ctg_depth_list

        line_num += 1

    #test_handle = open('/Users/songweizhi/Desktop/old.txt', 'w')

    bin_to_depth_dict = {}
    for genome_bin in sorted(bin_list):

        current_bin_ctg_list = bin_to_ctg_dict[genome_bin]
        current_bin_overall_depth = []
        n = 0
        while n < sample_num:

            current_bin_ctg_depth_weighted_sum = 0
            for current_bin_ctg in current_bin_ctg_list:
                current_bin_ctg_len = ctg_to_len_dict[current_bin_ctg]
                current_bin_ctg_depth = ctg_to_depth_dict[current_bin_ctg][n]
                current_bin_ctg_depth_weighted = current_bin_ctg_len*current_bin_ctg_depth
                current_bin_ctg_depth_weighted_sum += current_bin_ctg_depth_weighted

            current_bin_depth = current_bin_ctg_depth_weighted_sum/bin_to_size_dict[genome_bin]
            current_bin_overall_depth.append(float("{0:.2f}".format(current_bin_depth)))

            n += 1

        bin_to_depth_dict[genome_bin] = current_bin_overall_depth

  #      print('%s\t%s' % (genome_bin, current_bin_ctg_depth_weighted_sum))
 #       test_handle.write('%s\t%s\n' % (genome_bin, current_bin_ctg_depth_weighted_sum))
#    test_handle.close()

    # calculate bin depth
    bin_total_depth_by_sample = []
    n = 0
    while n < sample_num:

        current_sample_total_depth = 0
        for each_bin in bin_to_depth_dict:
            bin_depth = bin_to_depth_dict[each_bin]
            current_sample_total_depth += bin_depth[n]

        bin_total_depth_by_sample.append(float("{0:.3f}".format(current_sample_total_depth)))

        n += 1

    # get normalized bin depth
    bin_to_normalized_depth_dict = {}
    for genome_bin in sorted(bin_list):
        genome_bin_depth = bin_to_depth_dict[genome_bin]
        genome_bin_depth_normalized = []
        n = 0
        while n < len(genome_bin_depth):
            genome_bin_depth_normalized.append(float("{0:.2f}".format(genome_bin_depth[n]*100/bin_total_depth_by_sample[n])))
            n += 1

        bin_to_normalized_depth_dict[genome_bin] = genome_bin_depth_normalized
        #print('%s\t%s' % (genome_bin, genome_bin_depth_normalized))

    # write out abundance file
    out_put_txt_handle = open(output_abundance_txt, 'w')
    out_put_txt_handle.write('Bin_id\t%s\n' % '\t'.join(sample_name_list))
    for genome_bin in sorted(bin_list):
        genome_bin_depth_list = bin_to_normalized_depth_dict[genome_bin]
        genome_bin_depth_list_str = [str(i) for i in genome_bin_depth_list]
        out_put_txt_handle.write('%s\t%s\n' % (genome_bin, '\t'.join(genome_bin_depth_list_str)))
    out_put_txt_handle.close()

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

