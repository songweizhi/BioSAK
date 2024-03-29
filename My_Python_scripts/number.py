
##########################################################################################

num_1 = float("{0:.2f}".format(123.456789))  # 123.45
num_2 = float("{0:.3f}".format(123.456789))  # 123.456


# try to int, if not float
num_str = '123.456789'
try:
    pwy_cpl = int(num_str)
except:
    pwy_cpl = float(num_str)


def norm_num_list(num_list, norm_by, decimal_to_keep):
    num_list_after_norm = []
    for each_num in num_list:
        each_num_after_norm = float(("{0:.%sf}" % decimal_to_keep).format(each_num/norm_by))
        num_list_after_norm.append(each_num_after_norm)
    return num_list_after_norm


##########################################################################################

def numericable_str(my_string):

    is_numericable_str = False
    try:
        my_numeric_value = int(my_string)
        is_numericable_str = True
    except ValueError:
        try:
            my_numeric_value = float(my_string)
            is_numericable_str = True
        except ValueError:
            is_numericable_str = False

    return is_numericable_str


################################################### Distance matrix ####################################################

'''
Format_of_symmetric_distance_matrix_file (tab separated)

    A   B   C
A   1   2   3
B   4   5   6
C   7   8   9
'''

import pandas as pd

dm_file = 'demo_files/symmetric_distance_matrix.tab'

# read in symmetric distance matrix file
df = pd.read_table(dm_file, header=0, index_col=0)

# df to arrary
df_in_arrary = df.to_numpy()
# [[1 2 3]
#  [4 5 6]
#  [7 8 9]]

# get column name list
column_name_list = list(df.columns)
# ['A', 'B', 'C']

# get row name list
row_name_list = list(df.index)
# ['A', 'B', 'C']

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


################################### get best number of clusters ##################################

import kmedoids
import pandas as pd
from sklearn.metrics import silhouette_score


def get_num_of_clusters(dm_file):

    df = pd.read_table(dm_file, header=0, index_col=0)
    df_in_arrary = df.to_numpy()

    max_score = -1
    best_cluster_num = None
    for cluster_num in range(2, df.shape[0]):

        # perform clustering
        clusters = kmedoids.fasterpam(df, cluster_num)
        label_list_1_based = [int((i + 1)) for i in clusters.labels]

        # calculate silhouette score
        score = silhouette_score(df_in_arrary, label_list_1_based)
        if score > max_score:
            max_score = score
            best_cluster_num = cluster_num

    if max_score < 0.55:
        best_cluster_num = 1

    return best_cluster_num, max_score


dm_file = 'demo_files/edge_20_wd39_dm.tab'
cluster_num, score = get_num_of_clusters(dm_file)
print('Best\t%s\t%s' % (cluster_num, score))
# Best	2	0.8067079676284064


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


##########################################################################################

note = '''

# ln in numpy: np.log()
print(np.log(2.7182780911361077))
# 0.9999986251147811

'''
