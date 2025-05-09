import os
import argparse
from Bio import SeqIO
from datetime import datetime
from BioSAK.global_functions import time_format


select_seq_usage = '''
========================= select_seq example commands =========================

BioSAK select_seq -i input.fa -id seq_id.txt -o output.fa
BioSAK select_seq -i input.fa -id seq_id.txt -o output.fa -exclude

# ID file format: 
one id per line, do not include great than (>) symbol.

===============================================================================
'''

def select_seq(args):

    # read in argument
    seq_file        = args['i']
    id_file         = args['id']
    output_file     = args['o']
    exclude_seqs    = args['exclude']
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
    if exclude_seqs is False:
        print(datetime.now().strftime(time_format) + 'Extracting sequences in %s' % id_file)
    else:
        print(datetime.now().strftime(time_format) + 'Extracting sequences except those in %s' % id_file)

    # extract sequences
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        if exclude_seqs is False:
            if seq_id in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')

        else:
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

    select_seq_parser = argparse.ArgumentParser(usage=select_seq_usage)
    select_seq_parser.add_argument('-i',       required=True,                          help='fasta file')
    select_seq_parser.add_argument('-id',       required=True,                          help='sequence id,one id per line')
    select_seq_parser.add_argument('-o',       required=True,                          help='output file')
    select_seq_parser.add_argument('-exclude', required=False, action="store_true",    help='specify to extract sequences except those in -id')
    select_seq_parser.add_argument('-fq',      required=False, action="store_true",    help='in fastq format, default: fa')
    select_seq_parser.add_argument('-oneline', required=False, action="store_true",    help='put sequence in single line')
    args = vars(select_seq_parser.parse_args())
    select_seq(args)
