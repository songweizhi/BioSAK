################################## re-index column order #################################

import pandas
import numpy as np

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
