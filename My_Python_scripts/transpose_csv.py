import pandas as pd


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


file_in =  '/Users/songweizhi/Desktop/test.csv'
file_out = '/Users/songweizhi/Desktop/test_T.csv'
sep_symbol = '\t'    # ',' or '\t'
column_name_pos = 0  # set first row as column name
row_name_pos = 0     # set first column as row name
transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos)

