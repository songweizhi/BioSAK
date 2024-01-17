import argparse
from Bio import SeqIO

merge_seq_usage = '''
===================== merge_seq example commands =====================

BioSAK merge_seq -1 ctgs1.fa -2 ctgs2.fa -o combined_ctgs.fa -sl
BioSAK merge_seq -1 reads1.fq -2 reads2.fq -o combined_reads.fq -fq

======================================================================
'''


def merge_seq(args):

    # read in argument
    seq_file_1  = args['1']
    seq_file_2  = args['2']
    output_file = args['o']
    in_fastq    = args['fq']
    one_line    = args['sl']

    file_format = 'fasta'
    if in_fastq is True:
        file_format = 'fastq'

    wrote_seq_set = set()
    output_file_handle = open(output_file, 'w')

    # read in seq_file_1
    for each_seq in SeqIO.parse(seq_file_1, file_format):
        if each_seq.id not in wrote_seq_set:
            if in_fastq is False:
                if one_line is False:
                    SeqIO.write(each_seq, output_file_handle, 'fasta')
                else:
                    SeqIO.write(each_seq, output_file_handle, 'fasta-2line')
            else:
                SeqIO.write(each_seq, output_file_handle, 'fastq')
            wrote_seq_set.add(each_seq.id)

    # read in seq_file_2
    for each_seq in SeqIO.parse(seq_file_2, file_format):
        if each_seq.id not in wrote_seq_set:
            if in_fastq is False:
                if one_line is False:
                    SeqIO.write(each_seq, output_file_handle, 'fasta')
                else:
                    SeqIO.write(each_seq, output_file_handle, 'fasta-2line')
            else:
                SeqIO.write(each_seq, output_file_handle, 'fastq')
            wrote_seq_set.add(each_seq.id)

    output_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(usage=merge_seq_usage)
    parser.add_argument('-1',       required=True,                          help='input file 1')
    parser.add_argument('-2',       required=True,                          help='input file 2')
    parser.add_argument('-o',       required=True,                          help='output sequence file')
    parser.add_argument('-fq',      required=False, action="store_true",    help='in fastq format, default: fasta')
    parser.add_argument('-sl',      required=False, action="store_true",    help='write out sequence in single line format')
    args = vars(parser.parse_args())
    merge_seq(args)
