import argparse
from Bio import AlignIO


fa2phy_usage = '''
======= fa2phy example commands =======

TreeSAK fa2phy -i msa.fa -o msa.phy

=======================================
'''


def fa2phy(args):

    fasta_in  = args['i']
    phy_out   = args['o']

    alignment = AlignIO.read(fasta_in, 'fasta')

    max_seq_id_len = 0
    for each_seq in alignment:
        seq_id_len = len(each_seq.id)
        if seq_id_len > max_seq_id_len:
            max_seq_id_len = seq_id_len

    with open(phy_out, 'w') as msa_out_handle:
        msa_out_handle.write('%s %s\n' % (len(alignment), alignment.get_alignment_length()))
        for each_seq in alignment:
            seq_id = each_seq.id
            seq_id_with_space = '%s%s' % (seq_id, ' ' * (max_seq_id_len + 2 - len(seq_id)))
            msa_out_handle.write('%s%s\n' % (seq_id_with_space, str(each_seq.seq)))


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,   help='input MSA in fasta format')
    parser.add_argument('-o',      required=True,   help='output MSA in phylip format')
    args = vars(parser.parse_args())
    fa2phy(args)
