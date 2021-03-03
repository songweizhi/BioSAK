#!/usr/bin/env python
from __future__ import division
import os
import argparse
from Bio import SeqIO
from datetime import datetime
from BioSAK.global_functions import time_format
from BioSAK.global_functions import sep_path_basename_ext


select_seq_usage = '''
======================= select_seq example commands =======================

# Extract sequences with provided id
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.1.fa -option 1

# Extract sequences except those in seq_id.txt
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.0.fa -option 0

# seq_id.txt file format: one id per line, great than symbol exculded.

===========================================================================
'''

def select_seq(args):

    # read in argument
    seq_file        = args['seq']
    id_file         = args['id']
    select_option   = args['option']
    output_file     = args['out']

    # report
    if select_option == 1:
        print(datetime.now().strftime(time_format) + 'Extracting sequences in %s' % id_file)
    if select_option == 0:
        print(datetime.now().strftime(time_format) + 'Extracting sequences except those in %s' % id_file)

    # get provided id list
    seq_id_list = set()
    for seq_id in open(id_file):
        seq_id_list.add(seq_id.strip())

    # extract sequences
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        seq_id = seq_record.id

        if select_option == 1:
            if seq_id in seq_id_list:
                output_file_handle.write('>%s\n' % seq_record.id)
                output_file_handle.write('%s\n' % seq_record.seq)

        if select_option == 0:
            if seq_id not in seq_id_list:
                output_file_handle.write('>%s\n' % seq_record.id)
                output_file_handle.write('%s\n' % seq_record.seq)

    output_file_handle.close()

    print(datetime.now().strftime(time_format) + 'Extracted sequences exported to %s' % output_file)
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for select_seq
    parser.add_argument('-seq',       required=True,            help='sequence file')
    parser.add_argument('-id',        required=True,            help='sequence ids,one id per line')
    parser.add_argument('-option',    required=True, type=int,  help='choose from 0 and 1')
    parser.add_argument('-out',       required=True,            help='output file')

    args = vars(parser.parse_args())
    select_seq(args)
