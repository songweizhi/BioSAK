import os
import random
import argparse
from Bio import SeqIO


################################################# input #################################################

parser = argparse.ArgumentParser(description='', add_help=True)
required = parser.add_argument_group('required arguments')

required.add_argument('-i', dest='IN', nargs='?', required=True,  type=str, help='input reads file')
required.add_argument('-p', dest='Percent', nargs='?', required=True,  type=float, help='subsample percentage, e.g. 0.5 for 50%')
required.add_argument('-o', dest='OUT', nargs='?', required=True,  type=str, help='output reads file')

args = vars(parser.parse_args())
input_reads_file = args['IN']
sample_percent = args['Percent']
output_reads_file = args['OUT']

#########################################################################################################

# get total reads number
total_reads_num = 0
for each_read in SeqIO.parse(input_reads_file, 'fasta'):
    total_reads_num += 1


# get the number of reads to keep
reads_num_to_keep = round(total_reads_num * sample_percent)


# index all reads, start from 0
index_num_list = list(range(total_reads_num))


# randomly select reads too keep
to_keep_index = random.sample(index_num_list, reads_num_to_keep)
to_keep_index_sorted = sorted(to_keep_index)


# write out reads
output_file_handle = open(output_reads_file, 'w')
index = 0
for read in SeqIO.parse(input_reads_file, 'fasta'):
    if index in to_keep_index_sorted:
        output_file_handle.write('>%s\n' % read.id)
        output_file_handle.write('%s\n' % str(read.seq))
    index += 1
output_file_handle.close()

