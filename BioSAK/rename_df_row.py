import os
import argparse


rename_df_row_usage = '''
======================== rename_df_row example command ========================

BioSAK rename_df_row -i df_in.txt -r rename.txt -o df_out.txt
BioSAK rename_df_row -i df_in.txt -r rename.txt -o df_out.txt -skip_top_row

Note: dataframe need to be tab separated.

===============================================================================
'''


def rename_df_row(args):

    df_in           = args['i']
    rename_txt      = args['r']
    df_out          = args['o']
    skip_top_row    = args['skip_top_row']

    # check input file
    if rename_txt is not None:
        if os.path.isfile(df_in) is False:
            print('%s not found, program exited!' % df_in)
            exit()

    if rename_txt is not None:
        if os.path.isfile(rename_txt) is False:
            print('%s not found, program exited!' % rename_txt)
            exit()

    # get genome rename dict
    rename_dict = dict()
    for line in open(rename_txt):
        line_split = line.strip().split()
        rename_dict[line_split[0]] = line_split[1]

    # rename row header
    line_index = 1
    df_out_handle = open(df_out, 'w')
    for each_line in open(df_in):
        each_line_split = each_line.strip().split('\t')
        row_header      = each_line_split[0]
        row_header_new  = rename_dict.get(row_header, row_header)
        if line_index == 1:
            if skip_top_row is True:
                df_out_handle.write(each_line)
            else:
                df_out_handle.write(each_line.replace(row_header, row_header_new))
        else:
            df_out_handle.write(each_line.replace(row_header, row_header_new))
        line_index += 1
    df_out_handle.close()


if __name__ == '__main__':

    rename_df_row_parser = argparse.ArgumentParser(usage=rename_df_row_usage)
    rename_df_row_parser.add_argument('-i',             required=True,                        help='input dataframe')
    rename_df_row_parser.add_argument('-r',             required=True,                        help='rename file, tab separated')
    rename_df_row_parser.add_argument('-o',             required=True,                        help='output dataframe')
    rename_df_row_parser.add_argument('-skip_top_row',  required=False, action="store_true",  help='skip the first row')
    args = vars(rename_df_row_parser.parse_args())
    rename_df_row(args)
