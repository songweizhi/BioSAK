import os
import argparse
import numpy as np
import pandas as pd


subset_df_usage = '''
============================ subset_df example commands ============================

BioSAK subset_df -i demo_df.txt -r row_id.txt -o df_subset.txt -rm0col
BioSAK subset_df -i demo_df.txt -c col_id.txt -o df_subset.txt -e -rm0col
BioSAK subset_df -i demo_df.txt -r row_id.txt -c col_id.txt -o df_subset.txt -e
BioSAK subset_df -i demo_df.txt -r row_id.txt -c col_id.txt -o df_subset.txt -b -m

# Note
1. If the "-c" or "-r" file contains multiple columns (tab separated), the 1st 
   column will be used. 
2. If you want to skip the 1st row in the "-c" or "-r" file, provide "-skip1row".
3. When the -sr/-sc options are not specified,The order of rows/columns in the 
   output will match that of the -r/-c files, respectively.
4. File name for the -col_desc is fixed, arCOGdef.tab for arCOG, ko00001.keg for KEGG, 
   and CAZyDB.fam-activities.txt for CAZy.

====================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def subset_df(args):

    file_in             = args['i']
    file_out            = args['o']
    rows_to_keep_file   = args['r']
    cols_to_keep_file   = args['c']
    df_separator        = args['s']
    in_binary           = args['b']
    zero_as_minus_one   = args['m']
    skip1row            = args['skip1row']
    exclude_rows_cols   = args['e']
    rm_zero_col         = args['rm0col']
    sort_row            = args['sr']
    sort_col            = args['sc']
    col_desc_txt        = args['col_desc']
    row_desc_txt        = args['row_desc']
    replace_0_with_na   = args['na']

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

    # read in col_desc_txt
    col_desc_dict = dict()
    if col_desc_txt is not None:
        if os.path.isfile(col_desc_txt) is False:
            print('%s not found, program exited!' % col_desc_txt)
            exit()
        col_desc_txt_name, _, _, _ = sep_path_basename_ext(col_desc_txt)
        if col_desc_txt_name == 'arCOGdef.tab':
            for each_line in open(col_desc_txt, encoding="ISO-8859-1"):
                each_line_split = each_line.strip().split('\t')
                fun_id   = each_line_split[0]
                fun_cat  = each_line_split[1]
                fun_desc = each_line_split[3]
                col_desc_dict[fun_id] = '%s__%s__%s' % (fun_id, fun_cat, fun_desc)
        elif col_desc_txt_name == 'ko00001.keg':
            for each_line in open(col_desc_txt):
                if each_line[0] in ['A', 'B', 'C', 'D']:
                    each_line_split = each_line.strip().split(' ')
                    if each_line[0] == 'A':
                        current_A_id = each_line_split[0]
                        current_A_desc = ' '.join(each_line_split[1:])
                        col_desc_dict[current_A_id] = '%s__%s' % (current_A_id, current_A_desc)
                    elif each_line[0] == 'B':
                        if len(each_line_split) > 1:
                            current_B_id = each_line_split[2]
                            current_B_desc = ' '.join(each_line_split[3:])
                            col_desc_dict[current_B_id] = '%s__%s' % (current_B_id, current_B_desc)
                    elif each_line[0] == 'C':
                        current_C_id = each_line_split[4]
                        current_C_desc = ' '.join(each_line_split[5:])
                        col_desc_dict[current_C_id] = '%s__%s' % (current_C_id, current_C_desc)
                    elif each_line[0] == 'D':
                        current_D_id = each_line_split[6]
                        current_D_desc = ' '.join(each_line_split[7:])
                        col_desc_dict[current_D_id] = '%s__%s' % (current_D_id, current_D_desc)
        elif col_desc_txt_name == 'CAZyDB.fam-activities.txt':
            for each_fam in open(col_desc_txt_name):
                each_fam_split = each_fam.strip().split('	  ')
                if len(each_fam_split) == 2:
                    fam_id = each_fam_split[0]
                    fam_activities = each_fam_split[1]
                    col_desc_dict[fam_id] = fam_activities
                    col_desc_dict[fam_id] = '%s__%s' % (fam_id, fam_activities)
        else:
            for each_line in open(col_desc_txt):
                each_line_split = each_line.strip().split('\t')
                col_desc_dict[each_line_split[0]] = '%s__%s' % (each_line_split[0], each_line_split[1])

    # read in row_desc_txt
    row_desc_dict = dict()
    if row_desc_txt is not None:
        if os.path.isfile(row_desc_txt) is False:
            print('%s not found, program exited!' % row_desc_txt)
            exit()
        row_desc_txt_name, _, _, _ = sep_path_basename_ext(row_desc_txt)

    ###################################### get the id of rows and cols to subset #######################################

    # put all row and col headers in list
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_row_header_list = df.index.values.tolist()
    df_col_header_list = df.columns.values.tolist()

    # read in rows_to_keep_file
    rows_found_list   = []
    rows_missing_list = []
    if rows_to_keep_file is not None:
        if os.path.isfile(rows_to_keep_file) is False:
            print('%s not found, program exited!' % rows_to_keep_file)
            exit()
        else:
            file_row_index = 0
            for each_r in open(rows_to_keep_file):
                ignore_current_row = False
                if (skip1row is True) and (file_row_index == 0):
                    ignore_current_row = True
                if ignore_current_row is False:
                    row_id = each_r.strip().split()[0]
                    if row_id in df_row_header_list:
                        if row_id not in rows_found_list:
                            rows_found_list.append(row_id)
                    else:
                        if row_id not in rows_missing_list:
                            rows_missing_list.append(row_id)
                file_row_index += 1

    # read in cols_to_keep_file
    cols_found_list = []
    cols_missing_list = []
    if cols_to_keep_file is not None:
        if os.path.isfile(cols_to_keep_file) is False:
            print('%s not found, program exited!' % cols_to_keep_file)
            exit()
        else:
            file_row_index = 0
            for each_c in open(cols_to_keep_file):
                ignore_current_row = False
                if (skip1row is True) and (file_row_index == 0):
                    ignore_current_row = True
                if ignore_current_row is False:
                    col_id = each_c.strip().split()[0]
                    if col_id in df_col_header_list:
                        if col_id not in cols_found_list:
                            cols_found_list.append(col_id)
                    else:
                        if col_id not in cols_missing_list:
                            cols_missing_list.append(col_id)
                file_row_index += 1

    # report
    if len(rows_missing_list) > 0:
        print('The following rows are missing from the dataframe:\n%s'    % ','.join(sorted(list(rows_missing_list))))

    if len(cols_missing_list) > 0:
        print('The following columns are missing from the dataframe:\n%s' % ','.join(sorted(list(cols_missing_list))))

    ####################################################################################################################

    rows_to_keep_list = []
    cols_to_keep_list = []
    if exclude_rows_cols is False:
        rows_to_keep_list = rows_found_list
        cols_to_keep_list = cols_found_list
    else:
        for each_row in df_row_header_list:
            if each_row not in rows_found_list:
                rows_to_keep_list.append(each_row)
        for each_col in df_col_header_list:
            if each_col not in cols_found_list:
                cols_to_keep_list.append(each_col)

    # sort rows and cols
    if sort_row is True:
        rows_to_keep_list = sorted(rows_to_keep_list)
    if sort_col is True:
        cols_to_keep_list = sorted(cols_to_keep_list)

    if len(rows_to_keep_list) == 0:
        if len(cols_to_keep_list) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep_list]
    else:
        if len(cols_to_keep_list) == 0:
            subset_df = df.loc[rows_to_keep_list, :]
        else:
            subset_df = df.loc[rows_to_keep_list, cols_to_keep_list]

    if in_binary is True:
        subset_df[subset_df <= 0] = 0
        subset_df[subset_df > 0] = 1

    # turn 0 to -1
    if zero_as_minus_one is True:
        subset_df[subset_df == 0] = -1

    if rm_zero_col is True:
        subset_df = subset_df.loc[:, (subset_df != 0).any(axis=0)]

    # rename columns
    if len(col_desc_dict) > 0:
        current_columns = subset_df.columns.tolist()
        new_columns = [col_desc_dict[i] for i in current_columns]
        subset_df.columns = new_columns

    # rename rows
    if len(row_desc_dict) > 0:
        pass
        # to be added

    if replace_0_with_na is True:
        subset_df.replace(0, 'NA', inplace=True)

    subset_df.to_csv(file_out, sep=sep_symbol)


if __name__ == '__main__':

    subset_df_parser = argparse.ArgumentParser(usage=subset_df_usage)
    subset_df_parser.add_argument('-i',         required=True,                          help='input file')
    subset_df_parser.add_argument('-o',         required=True,                          help='output file')
    subset_df_parser.add_argument('-c',         required=False, default=None,           help='columns to keep')
    subset_df_parser.add_argument('-r',         required=False, default=None,           help='rows to keep')
    subset_df_parser.add_argument('-s',         required=False, default='tab',          help='column separator, choose from tab and comma, default: tab')
    subset_df_parser.add_argument('-b',         required=False, action='store_true',    help='write out dataframe in 0/1 format')
    subset_df_parser.add_argument('-e',         required=False, action='store_true',    help='subset df by excluding the specified rows/columns')
    subset_df_parser.add_argument('-m',         required=False, action='store_true',    help='convert 0 to -1')
    subset_df_parser.add_argument('-rm0col',    required=False, action='store_true',    help='remove columns contain only 0')
    subset_df_parser.add_argument('-skip1row',  required=False, action='store_true',    help='skip the 1st row of the -c/-r file')
    subset_df_parser.add_argument('-sr',        required=False, action='store_true',    help='sort rows')
    subset_df_parser.add_argument('-sc',        required=False, action='store_true',    help='sort columns')
    subset_df_parser.add_argument('-col_desc',  required=False, default=None,           help='column description')
    subset_df_parser.add_argument('-row_desc',  required=False, default=None,           help='row description')
    subset_df_parser.add_argument('-na',        required=False, action='store_true',    help='replace 0 with NA')
    args = vars(subset_df_parser.parse_args())
    subset_df(args)
