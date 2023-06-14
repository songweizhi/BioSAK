import os
import argparse


metabat2maxbin_usage = '''
=============== metabat2maxbin example commands ===============

BioSAK metabat2maxbin -p maxbin2_depth -i metabat_depth.txt

===============================================================
'''

def metabat2maxbin(args):

    maxbin_depth_prefix = args['p']
    metabat_depth_txt   = args['i']

    to_write_dict = dict()
    line_index = 1
    for each_line in open(metabat_depth_txt):
        each_line_split = each_line.strip().split('\t')

        col_index = 0
        sample_index = 1
        for each_col in each_line_split:
            if (col_index >= 3) and ((col_index % 2) != 0):

                op_txt = '%s_%s.txt' % (maxbin_depth_prefix, sample_index)
                if op_txt not in to_write_dict:
                    to_write_dict[op_txt] = []

                if line_index >= 2:  # to ignore the header line
                    to_write_dict[op_txt].append('%s\t%s' % (each_line_split[0], each_col))

                sample_index += 1
            col_index += 1
        line_index += 1

    for each_file in to_write_dict:
        file_handle = open(each_file, 'w')
        file_handle.write('\n'.join(to_write_dict[each_file]))
        file_handle.close()


if __name__ == '__main__':

    metabat2maxbin_parser = argparse.ArgumentParser(usage=metabat2maxbin_usage)
    metabat2maxbin_parser.add_argument('-p',   required=True,  help='output prefix for maxbin depth file')
    metabat2maxbin_parser.add_argument('-i',   required=True,  help='metabat depth')
    args = vars(metabat2maxbin_parser.parse_args())
    metabat2maxbin(args)
