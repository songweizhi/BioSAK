############################################## subset_df ##############################################

import pandas as pd

def subset_df(file_in, rows_to_keep, cols_to_keep, sep_symbol, row_name_pos, column_name_pos, file_out):

    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)

    if len(rows_to_keep) == 0:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep]
    else:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[rows_to_keep, :]
        else:
            subset_df = df.loc[rows_to_keep, cols_to_keep]

    subset_df.to_csv(file_out, sep=sep_symbol)

rows_to_keep    = set()
cols_to_keep    = set()
sep_symbol      = '\t'
row_name_pos    = 0
column_name_pos = 0
