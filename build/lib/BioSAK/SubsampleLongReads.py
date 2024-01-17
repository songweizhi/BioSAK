import os
import random
import argparse
from Bio import SeqIO


SubsampleLongReads_usage = '''
====================== SubsampleLongReads example commands ======================

BioSAK SubsampleLongReads -i LongReads.fastq -s 5 -o op_dir_5 -fq
BioSAK SubsampleLongReads -i LongReads.fastq -s 1,5,10 -o op_dir_1-10 -fq
BioSAK SubsampleLongReads -i LongReads.fastq -s 1,25,50,75 -o op_dir_1-75 -fq

=================================================================================
'''


def SubsampleLongReads(args):

    seq_file        = args['i']
    subsample_str   = args['s']
    op_dir          = args['o']
    one_line        = args['oneline']
    in_fastq        = args['fq']

    # check if output dir exist
    if os.path.isdir(op_dir):
        print('The specified output directory already exist, program exited!')
        print(op_dir)
        exit()

    os.system('mkdir %s' % op_dir)

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    # read in read id
    seq_id_set = set()
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        seq_id_set.add(seq_id)

    subsample_dict = dict()
    for subsample_pct in subsample_str.split(','):
        read_num_to_subset = round(len(seq_id_set)*float(subsample_pct)/100)
        reads_to_extract = random.sample([i for i in seq_id_set], read_num_to_subset)
        for each_read in reads_to_extract:
            if each_read not in subsample_dict:
                subsample_dict[each_read] = {subsample_pct}
            else:
                subsample_dict[each_read].add(subsample_pct)

    # subsampling
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        if seq_id in subsample_dict:
            subsample_pct_set = subsample_dict[seq_id]
            for each_subsample_pct in subsample_pct_set:
                current_op_file = '%s/%s.%s' % (op_dir, each_subsample_pct, seq_in_format)
                with open (current_op_file, 'a') as current_op_file_handle:
                    if in_fastq is False:
                        if one_line is False:
                            SeqIO.write(seq_record, current_op_file_handle, 'fasta')
                        else:
                            SeqIO.write(seq_record, current_op_file_handle, 'fasta-2line')
                    else:
                        SeqIO.write(seq_record, current_op_file_handle, 'fastq')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=SubsampleLongReads_usage)
    parser.add_argument('-i',       required=True,                       help='input sequence file')
    parser.add_argument('-s',       required=True,                       help='separate subsample rates (1-100) with comma, e.g. 1,5,10')
    parser.add_argument('-o',       required=True,                       help='output directory')
    parser.add_argument('-fq',      required=False, action="store_true", help='in fastq format, default: fa')
    parser.add_argument('-oneline', required=False, action="store_true", help='put sequence in single line')
    args = vars(parser.parse_args())
    SubsampleLongReads(args)
