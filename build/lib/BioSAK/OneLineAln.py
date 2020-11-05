import argparse
from Bio import SeqIO


OneLineAln_usage = '''
========================= OneLineAln example commands =========================

BioSAK OneLineAln -in MarkerGenes.aln -out MarkerGenes_OneLine.aln
BioSAK OneLineAln -in MarkerGenes.aln -out MarkerGenes_OneLine.aln -upper

===============================================================================
'''

def OneLineAln(args):

    aln_in_fasta     = args['in']
    aln_out_one_line = args['out']
    to_uppercase     = args['upper']

    # get longest_seq_id
    longest_seq_id = 0
    for seq in SeqIO.parse(aln_in_fasta, 'fasta'):
        if len(seq.id) > longest_seq_id:
            longest_seq_id = len(seq.id)

    # write out in new format
    aln_in_one_line_handle = open(aln_out_one_line, 'w')
    for seq in SeqIO.parse(aln_in_fasta, 'fasta'):
        seq_id_polished = seq.id + (longest_seq_id - len(seq.id))*' '
        if to_uppercase is True:
            aln_in_one_line_handle.write('%s\t%s\n' % (seq_id_polished, str(seq.seq).upper()))
        else:
            aln_in_one_line_handle.write('%s\t%s\n' % (seq_id_polished, str(seq.seq)))

    aln_in_one_line_handle.close()


if __name__ == '__main__':

    OneLineAln_parser = argparse.ArgumentParser()

    # arguments for rename_seq_parser
    OneLineAln_parser.add_argument('-in',     required=True,                        help='input MSA in fasta format')
    OneLineAln_parser.add_argument('-out',    required=False, default=None,         help='output file')
    OneLineAln_parser.add_argument('-upper',  required=False, action='store_true',  help='turn to uppercase')

    args = vars(OneLineAln_parser.parse_args())

    OneLineAln(args)

