import os
import argparse
from Bio import SeqIO


fq2fa_usage = '''
========== fq2fa example commands ==========

BioSAK fq2fa -i reads_R1.fq -o reads_R1.fa
BioSAK fq2fa -i reads_R2.fq -o reads_R2.fa

============================================
'''

def fq2fa(args):

    fq_in   = args['i']
    fa_out  = args['o']

    fa_out_handle = open(fa_out, 'w')
    for each_long_read in SeqIO.parse(fq_in, 'fastq'):
        fa_out_handle.write('>%s\n%s\n' % (each_long_read.id, each_long_read.seq))
    fa_out_handle.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='input fastq')
    parser.add_argument('-o', required=True, help='output fasta')
    args = vars(parser.parse_args())
    fq2fa(args)
