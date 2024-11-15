#!/usr/bin/python3
import argparse
from Bio import SeqIO


def export_seq(seq_record, single_line, output_handle):
    if single_line == 0:
        SeqIO.write(seq_record, output_handle, 'fasta')
    if single_line == 1:
        output_handle.write('>%s\n%s\n' % (seq_record.id, str(seq_record.seq)))


parser = argparse.ArgumentParser()

parser.add_argument('-in',
                    required=True,
                    help='input file')

parser.add_argument('-out',
                    required=True,
                    help='output file')

parser.add_argument('-step',
                    required=True,
                    type=int,
                    help='step')

parser.add_argument('-single_line',
                    action="store_true",
                    required=False,
                    help='to export sequences in single line')

args = vars(parser.parse_args())
input_file = args['in']
output_file = args['out']
step = args['step']
single_line = args['single_line']


n = 0
m = 0
output_handle = open(output_file, 'w')
for each_seq_record in SeqIO.parse(input_file, 'fasta'):
    n += 1
    if (n == 2*(step*m) + 1) or (n == 2*(step*m) + 2):
        export_seq(each_seq_record, single_line, output_handle)
        if n == 2*(step*m) + 2:
            m += 1
output_handle.close()
