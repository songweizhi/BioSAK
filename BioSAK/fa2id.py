import os
import argparse
from Bio import SeqIO

fa2id_usage = '''
============== fa2id example command ==============

BioSAK fa2id -i gnm.fa -o ctg_id.txt
BioSAK fa2id -i reads.fastq -o reads_id.txt -fq

===================================================
'''


def fa2id(args):

    fasta_in   = args['i']
    output_txt = args['o']
    in_fastq   = args['fq']

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    output_txt_handle = open(output_txt, 'w')
    for each_seq in SeqIO.parse(fasta_in, seq_in_format):
        output_txt_handle.write(each_seq.id + '\n')
    output_txt_handle.close()


if __name__ == '__main__':

    fa2id_usage_parser = argparse.ArgumentParser(usage=fa2id_usage)
    fa2id_usage_parser.add_argument('-i',  required=True, help='input fasta/fastq file')
    fa2id_usage_parser.add_argument('-o',  required=True, help='output txt')
    fa2id_usage_parser.add_argument('-fq', required=False, action="store_true",  help='in fastq format, default: fasta')
    args = vars(fa2id_usage_parser.parse_args())
    fa2id(args)
