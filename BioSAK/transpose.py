import argparse
import pandas as pd


transpose_usage = '''
============ transpose example commands ============

BioSAK transpose -i demo_df.txt -o demo_df_T.txt

# Input data matrix need to be tab separated

====================================================
'''


def transpose(args):

    file_in  = args['i']
    file_out = args['o']

    df = pd.read_csv(file_in, sep='\t', header=0, index_col=0)
    df_t = df.T
    df_t.to_csv(file_out, sep='\t')


if __name__ == '__main__':

    transpose_parser = argparse.ArgumentParser(usage=transpose_usage)
    transpose_parser.add_argument('-i',     required=True,                          help='input file')
    transpose_parser.add_argument('-o',     required=True,                          help='output file')
    args = vars(transpose_parser.parse_args())
    transpose(args)
