#!/usr/bin/env python3

import os
import argparse
from BioSAK.global_functions import sep_path_basename_ext


fa2tree_usage = '''
============ fa2tree example commands ============

module load mafft/7.407
module load fasttree/2.1.11
BioSAK fa2tree -seq 16S.fa -t 12

==================================================
'''


def fa2tree(args):

    seq_file    = args['seq']
    thread_num  = args['t']

    seq_file_path, seq_file_basename, seq_file_ext = sep_path_basename_ext(seq_file)

    aln_file  = '%s/%s.aln'     % (seq_file_path, seq_file_basename)
    tree_file = '%s/%s.newick'  % (seq_file_path, seq_file_basename)

    cmd_mafft    = 'mafft --quiet --maxiterate 1000 --adjustdirection --thread %s --globalpair %s > %s' % (thread_num, seq_file, aln_file)
    print(cmd_mafft)
    os.system(cmd_mafft)

    cmd_fasttree = 'fasttree -nt -quiet %s > %s' % (aln_file, tree_file)
    print(cmd_fasttree)
    os.system(cmd_fasttree)

    print('Tree exported to %s' % tree_file)


if __name__ == '__main__':

    fa2tree_parser = argparse.ArgumentParser(usage=fa2tree_usage)
    fa2tree_parser.add_argument('-seq',  required=True,                        help='sequence file')
    fa2tree_parser.add_argument('-t',    required=False, type=int, default=1,  help='number of threads, default: 1')
    args = vars(fa2tree_parser.parse_args())
    fa2tree(args)
