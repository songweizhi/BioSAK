import os
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help='sequence file')
args = vars(parser.parse_args())
file_in     = args['in']


file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
file_out = '%s/%s_stats.txt' % (file_in_path, file_in_basename)

total_number = 0
total_length = 0
for seq_record in SeqIO.parse(file_in, 'fastq'):
    total_length += len(seq_record.seq)
    total_number += 1


file_out_handle = open(file_out, 'w')
file_out_handle.write('%s\t%s\t%sbp\n' % (file_in, total_number, total_length))
file_out_handle.close()
