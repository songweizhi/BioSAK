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


# read in argument
parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help='input fastq file')
args = vars(parser.parse_args())
fq_in = args['in']


# define output fasta file name
fq_in_path, fq_in_basename, fq_in_extension = sep_path_basename_ext(fq_in)
fa_out = '%s/%s.fa' % (fq_in_path, fq_in_basename)


# extract sequences
fa_out_handle = open(fa_out, 'w')
for seq_record in SeqIO.parse(fq_in, 'fastq'):
    seq_record_sequence = str(seq_record.seq)
    seq_record_description = seq_record.description
    fa_out_handle.write('>%s\n' % seq_record_description)
    fa_out_handle.write('%s\n' % seq_record_sequence)
fa_out_handle.close()
