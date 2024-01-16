import os
import argparse
from Bio import SeqIO
from datetime import datetime
from BioSAK.global_functions import time_format


select_seq_usage = '''
======================= select_seq example commands =======================

# Extract sequences with provided id
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.1.fa -option 1 
BioSAK select_seq -seq reads.fastq -id read_id.txt -out output.1.fq -option 1 -fq

# Extract sequences except those in seq_id.txt
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.0.fa -option 0
BioSAK select_seq -seq ctg.fasta -id seq_id.txt -out output.0.fa -option 0 -oneline

# seq_id.txt file format: one id per line, great than symbol excluded.

===========================================================================
'''

def select_seq(args):

    # read in argument
    seq_file        = args['seq']
    id_file         = args['id']
    select_option   = args['option']
    output_file     = args['out']
    one_line        = args['oneline']
    in_fastq        = args['fq']

    # get provided id list
    print(datetime.now().strftime(time_format) + 'Read in sequence id in %s' % id_file)
    seq_id_list = set()
    for seq_id in open(id_file):
        seq_id_list.add(seq_id.strip())

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    # report
    if select_option == 1:
        print(datetime.now().strftime(time_format) + 'Extracting sequences in %s' % id_file)
    if select_option == 0:
        print(datetime.now().strftime(time_format) + 'Extracting sequences except those in %s' % id_file)

    # extract sequences
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        if select_option == 1:
            if seq_id in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')

        if select_option == 0:
            if seq_id not in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')

    output_file_handle.close()

    print(datetime.now().strftime(time_format) + 'Sequences exported to %s' % output_file)
    print(datetime.now().strftime(time_format) + 'Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=select_seq_usage)
    parser.add_argument('-seq',       required=True,                        help='sequence file')
    parser.add_argument('-id',        required=True,                        help='sequence ids,one id per line')
    parser.add_argument('-option',    required=True, type=int,              help='choose from 0 and 1')
    parser.add_argument('-out',       required=True,                        help='output file')
    parser.add_argument('-fq',        required=False, action="store_true",  help='in fastq format, default: fa')
    parser.add_argument('-oneline',   required=False, action="store_true",  help='put sequence in single line')
    args = vars(parser.parse_args())
    select_seq(args)
