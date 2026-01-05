import argparse
from Bio import SeqIO


rm_short_usage = '''
========= rm_short example commands =========

BioSAK rm_short -f -i in.fa -o out.fa -l 50

=============================================
'''


def rm_short(args):

    fa_in   = args['i']
    fa_out  = args['o']
    min_len = args['l']

    op_file_handle = open(fa_out, 'w')
    for each_seq in SeqIO.parse(fa_in, 'fasta'):
        if len(each_seq.seq) >= min_len:
            op_file_handle.write('>%s\n' % each_seq.id)
            op_file_handle.write('%s\n' % each_seq.seq)
    op_file_handle.close()


if __name__ == '__main__':

    rm_short_parser = argparse.ArgumentParser(usage=rm_short_usage)
    rm_short_parser.add_argument('-i',   required=True,             help='input fasta file')
    rm_short_parser.add_argument('-l',   required=True,type=int,    help='length cutoff in bp')
    rm_short_parser.add_argument('-o',   required=True,             help='output fasta file')
    args = vars(rm_short_parser.parse_args())
    rm_short(args)
