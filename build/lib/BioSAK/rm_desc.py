import argparse
from Bio import SeqIO


rm_desc_usage = '''
====== rm_desc example commands ======

BioSAK rm_desc -f -i in.fa -o out.fa

======================================
'''


def rm_desc(args):

    fa_in   = args['i']
    fa_out  = args['o']

    op_file_handle = open(fa_out, 'w')
    for each_seq in SeqIO.parse(fa_in, 'fasta'):
        op_file_handle.write('>%s\n' % each_seq.id)
        op_file_handle.write('%s\n' % each_seq.seq)
    op_file_handle.close()


if __name__ == '__main__':

    rm_desc_parser = argparse.ArgumentParser(usage=rm_desc_usage)
    rm_desc_parser.add_argument('-i',   required=True,  help='input fasta file')
    rm_desc_parser.add_argument('-o',   required=True,  help='output fasta file')
    args = vars(rm_desc_parser.parse_args())
    rm_desc(args)
