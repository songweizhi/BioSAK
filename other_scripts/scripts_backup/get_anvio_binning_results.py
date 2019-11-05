import os
import glob
import argparse
from Bio import SeqIO


##################################################### CONFIGURATION ####################################################

parser = argparse.ArgumentParser()

parser.add_argument('-f',
                    required=True,
                    help='path to bin folder')

parser.add_argument('-x',
                    required=True,
                    help='bin file extension')

args = vars(parser.parse_args())
bin_folder = args['f']
bin_ext = args['x']

if bin_folder[-1] == '/':
    bin_folder = bin_folder[:-1]

########################################################################################################################

output_file = 'binning_results.txt'

bin_re = '%s/*.%s' % (bin_folder, bin_ext)
bin_file_list = [os.path.basename(file_name) for file_name in glob.glob(bin_re)]


bin_folder_p = '/'.join(bin_folder.split('/')[:-1])
if bin_folder_p == '':
    bin_folder_p = '.'
pwd_output_file = '%s/%s' % (bin_folder_p, output_file)


output_file_handle = open(pwd_output_file, 'w')
for each_bin in bin_file_list:
    each_bin_basename, ext = os.path.splitext(each_bin)
    pwd_each_bin = '%s/%s' % (bin_folder, each_bin)
    for each_ctg in SeqIO.parse(pwd_each_bin, 'fasta'):
        for_write = '%s\t%s\n' % (each_ctg.id, each_bin_basename)
        output_file_handle.write(for_write)

output_file_handle.close()

