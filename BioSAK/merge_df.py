import argparse
import pandas as pd


merge_df_usage = '''
================ merge_df example commands ================

BioSAK merge_df -1 df_1.txt -2 df_2.txt -o merged_df.txt

/Library/Frameworks/Python.framework/Versions/3.10/bin/python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/merge_df.py -1 df_1.txt -2 df_2.txt -o merged_df.txt

===========================================================
'''


def merge_df(args):

    df1_txt             = args['1']
    df2_txt             = args['2']
    file_out            = args['o']
    df_separator        = args['s']
    column_name_pos     = 0
    row_name_pos        = 0

    # setup separator
    if df_separator in ['tab', 'Tab', 'TAB']:
        sep_symbol = '\t'
    elif df_separator in ['comma', 'Comma', 'COMMA']:
        sep_symbol = ','
    else:
        print('Please specify separator as either tab or comma, program exited!')
        exit()

    df1 = pd.read_csv(df1_txt, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df2 = pd.read_csv(df2_txt, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)

    frames = [df1, df2]
    merged_df = pd.concat(frames, axis=1)

    merged_df.to_csv(file_out, sep=sep_symbol, na_rep='na')


if __name__ == '__main__':

    merge_df_parser = argparse.ArgumentParser(usage=merge_df_usage)
    merge_df_parser.add_argument('-1',     required=True,                          help='the first input dataframe')
    merge_df_parser.add_argument('-2',     required=True,                          help='the second input dataframe')
    merge_df_parser.add_argument('-o',     required=True,                          help='output file')
    merge_df_parser.add_argument('-s',     required=False, default='tab',          help='column separator, choose from tab and comma, default: tab')
    args = vars(merge_df_parser.parse_args())
    merge_df(args)
