import os
import argparse
import pandas as pd


subset_df_usage = '''
================================== subset_df example commands ==================================

BioSAK subset_df -i demo_df.txt -r rows_to_keep.txt -o df_subset.txt
BioSAK subset_df -i demo_df.txt -c cols_to_keep.txt -o df_subset.txt
BioSAK subset_df -i demo_df.txt -r rows_to_keep.txt -c cols_to_keep.txt -o df_subset.txt
BioSAK subset_df -i demo_df.txt -r rows_to_keep.txt -c cols_to_keep.txt -o df_subset.txt -b -m

================================================================================================
'''


def subset_df(args):

    file_in             = args['i']
    file_out            = args['o']
    rows_to_keep_file   = args['r']
    cols_to_keep_file   = args['c']
    df_separator        = args['s']
    in_binary           = args['b']
    zero_as_minus_one   = args['m']
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

    ###################################### get the id of rows and cols to subset #######################################

    # put all row and col headers in list
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_row_header_list = df.index.values.tolist()
    df_col_header_list = df.columns.values.tolist()

    # read in rows_to_keep_file
    rows_to_keep_set = set()
    rows_missing_set = set()
    if os.path.isfile(rows_to_keep_file) is True:
        for each_r in open(rows_to_keep_file):
            row_id = each_r.strip().split()[0]
            if row_id in df_row_header_list:
                rows_to_keep_set.add(row_id)
            else:
                rows_missing_set.add(row_id)

    # read in cols_to_keep_file
    cols_to_keep_set = set()
    cols_missing_set = set()
    if os.path.isfile(cols_to_keep_file) is True:
        for each_c in open(cols_to_keep_file):
            col_id = each_c.strip().split()[0]
            if col_id in df_col_header_list:
                cols_to_keep_set.add(col_id)
            else:
                cols_missing_set.add(col_id)

    # report
    if len(rows_missing_set) > 0:
        print('The following rows are missing from the dataframe:\n%s'    % ','.join(sorted(list(rows_missing_set))))

    if len(cols_missing_set) > 0:
        print('The following columns are missing from the dataframe:\n%s' % ','.join(sorted(list(cols_missing_set))))

    ####################################################################################################################

    # turn sets into lists
    rows_to_keep_list_sorted = sorted(list(rows_to_keep_set))
    cols_to_keep_list_sorted = sorted(list(cols_to_keep_set))

    if len(rows_to_keep_list_sorted) == 0:
        if len(cols_to_keep_list_sorted) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep_list_sorted]
    else:
        if len(cols_to_keep_list_sorted) == 0:
            subset_df = df.loc[rows_to_keep_list_sorted, :]
        else:
            subset_df = df.loc[rows_to_keep_list_sorted, cols_to_keep_list_sorted]

    if in_binary is True:
        subset_df[subset_df <= 0] = 0
        subset_df[subset_df > 0] = 1

    # turn 0 to -1
    if zero_as_minus_one is True:
        subset_df[subset_df == 0] = -1

    subset_df.to_csv(file_out, sep=sep_symbol)


if __name__ == '__main__':

    subset_df_parser = argparse.ArgumentParser(usage=subset_df_usage)
    subset_df_parser.add_argument('-i',     required=True,                          help='input file')
    subset_df_parser.add_argument('-o',     required=True,                          help='output file')
    subset_df_parser.add_argument('-r',     required=False, default='',             help='header of rows to keep')
    subset_df_parser.add_argument('-c',     required=False, default='',             help='header of columns to keep')
    subset_df_parser.add_argument('-s',     required=False, default='tab',          help='column separator, choose from tab and comma, default: tab')
    subset_df_parser.add_argument('-b',     required=False, action='store_true',    help='write out dataframe in 0/1 format')
    subset_df_parser.add_argument('-m',     required=False, action='store_true',    help='convert 0 to -1')
    args = vars(subset_df_parser.parse_args())
    subset_df(args)
