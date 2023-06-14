import os
import argparse


metabat2concoct_usage = '''
================= metabat2concoct example commands =================

BioSAK metabat2concoct -i metabat_depth.txt -o concoct_depth.txt

====================================================================
'''

def metabat2concoct(args):

    metabat_depth_txt       = args['i']
    concoct_depth_txt       = args['o']

    concoct_depth_txt_handle = open(concoct_depth_txt, 'w')
    for each_line in open(metabat_depth_txt):
        each_line_split = each_line.strip().split('\t')

        col_index = 0
        cols_to_keep = []
        for each_col in each_line_split:

            if col_index == 0:
                cols_to_keep.append(each_col)
            elif col_index >= 3:
                if (col_index % 2) != 0:
                    cols_to_keep.append(each_col)
            col_index += 1
        concoct_depth_txt_handle.write('\t'.join(cols_to_keep) + '\n')
    concoct_depth_txt_handle.close()


if __name__ == '__main__':

    metabat2concoct_parser = argparse.ArgumentParser(usage=metabat2concoct_usage)
    metabat2concoct_parser.add_argument('-i',   required=True,  help='metabat depth')
    metabat2concoct_parser.add_argument('-o',   required=True,  help='concoct depth')
    args = vars(metabat2concoct_parser.parse_args())
    metabat2concoct(args)
