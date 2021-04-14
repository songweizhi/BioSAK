
note = '''

# ln in numpy: np.log()
print(np.log(2.7182780911361077))
# 0.9999986251147811

'''


##########################################################################################

num = 123.456789
num_1 = float("{0:.2f}".format(num))
num_2 = float("{0:.3f}".format(num))
print(num_1)    # 123.46
print(num_2)    # 123.457


def norm_num_list(num_list, norm_by, decimal_to_keep):

    num_list_after_norm = []
    for each_num in num_list:
        each_num_after_norm = float(("{0:.%sf}" % decimal_to_keep).format(each_num/norm_by))
        num_list_after_norm.append(each_num_after_norm)

    return num_list_after_norm


################################### cluster number list ##################################

def cluster(data, maxgap, min_size):

    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])

    groups_qualified = []
    for each_cluster in groups:
        if len(each_cluster) >= min_size:
            groups_qualified.append(each_cluster)

    return groups_qualified


######################### get_malthusian_selection_coefficients ##########################

# Equation for calculating mij in reference:
# Fitness effects of advantageous mutations in evolving Escherichia coli populations

import numpy as np


def get_malthusian_selection_coefficients(mutant_freq_start, mutant_freq_end, generation_num):

    malthusian_selection_coefficients = np.log((mutant_freq_end / mutant_freq_start) / ((1 - mutant_freq_end) / (1 - mutant_freq_start))) / generation_num

    return malthusian_selection_coefficients


mutant_freq_start   = 0.2
mutant_freq_end     = 0.8
generation_num      = 50
Malthusian_Selection_Coefficients = get_malthusian_selection_coefficients(mutant_freq_start, mutant_freq_end, generation_num)
print(Malthusian_Selection_Coefficients)


##################################### get percentile #####################################

import numpy as np

num_list = [1,2,3,4,5,6,7,8,9,10]
num_array = np.array(num_list)

percentile_25 = np.percentile(num_array, 25)
percentile_75 = np.percentile(num_array, 75)

print(percentile_25)	# 3.25
print(percentile_75)	# 7.75


###################################### transpose_csv #####################################

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



################################## re-index column order #################################

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


##########################################################################################

import os
import numpy as np
import pandas as pd
from scipy import stats
import Kelp_HGT_config


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


##########################################################################################

import math

print(int(math.ceil(11/float(3))))    # 4

# The divmod() returns a pair of numbers (a tuple) consisting of quotient q and remainder r
print(divmod(20, 3))    # (6, 2)


import numpy as np
from scipy import stats

a = [1, 2, 3, 3, 3, 3, 4, 5]
b = [5, 6, 7, 4, 6, 7, 4, 6, 6, 1]

print(np.sum(a))
print(np.mean(a))
print(np.std(a))
print('\n')
print(np.sum(b))
print(np.mean(b))
print(np.std(b))
print('\n')
p = stats.ttest_ind(a,b)
print(p)



