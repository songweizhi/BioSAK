
########################################################################################################################

num = 123.456789
num_1 = float("{0:.2f}".format(num))
num_2 = float("{0:.3f}".format(num))

print(num_1)    # 123.46
print(num_2)    # 123.457


#################################################### get percentile ####################################################

import numpy as np

num_list = [1,2,3,4,5,6,7,8,9,10]
num_array = np.array(num_list)

percentile_25 = np.percentile(num_array, 25)
percentile_75 = np.percentile(num_array, 75)

print(percentile_25)	# 3.25
print(percentile_75)	# 7.75


#################################################### transpose_csv #####################################################

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



################################################# re-index column order ################################################

import pandas

# turn list into arrary
category_num_arrary = np.array([[1,2,3], [4,5,6], [7,8,9]])

# add row and column name to dataframe
category_num_df = pandas.DataFrame(category_num_arrary, index=['Sample_1', 'Sample_2', 'Sample_3'], columns=['B', 'A', 'C'])
print(category_num_df)
# re-index column order
category_num_df = category_num_df.reindex(['A', 'B', 'C'], axis=1)
print(category_num_df)

# write out
#category_num_df.to_csv('Output.csv')

'''
          B  A  C                    A  B  C
Sample_1  1  2  3          Sample_1  2  1  3
Sample_2  4  5  6    -->   Sample_2  5  4  6
Sample_3  7  8  9          Sample_3  8  7  9

'''


########################################################################################################################

import math

print(int(math.ceil(11/float(3))))    # 4

# The divmod() returns a pair of numbers (a tuple) consisting of quotient q and remainder r
print(divmod(20, 3))    # (6, 2)
