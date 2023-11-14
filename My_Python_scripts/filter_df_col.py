##################################### filter dataframe columns #####################################

import pandas as pd


def filter_df_col(df_in, min_value, min_value_min_num):

    df_in_copy = df_in.copy(deep=True)

    for col in df_in_copy.columns:
        count = (df_in_copy[col] >= min_value).sum()
        if count < min_value_min_num:
            print(count)
            df_in_copy.drop(col, axis=1, inplace=True)

    return df_in_copy


df = pd.DataFrame({
    'col1': [1, 1, 1, 1, 1, 1, 0],
    'col2': [1, 1, 1, 1, 1, 0, 0],
    'col3': [1, 1, 1, 1, 0, 0, 0],
    'col4': [1, 1, 1, 0, 0, 0, 0],
    'col5': [1, 1, 0, 0, 0, 0, 0],
    'col6': [1, 0, 0, 0, 0, 0, 0],
    'col7': [0, 0, 0, 0, 0, 0, 0]
})

min_value           = 1
min_value_min_num   = 5
df_filtered         = filter_df_col(df, min_value, min_value_min_num)


'''
print(df)

   col1  col2  col3  col4  col5  col6  col7
0     1     1     1     1     1     1     0
1     1     1     1     1     1     0     0
2     1     1     1     1     0     0     0
3     1     1     1     0     0     0     0
4     1     1     0     0     0     0     0
5     1     0     0     0     0     0     0
6     0     0     0     0     0     0     0

print(df_filtered)

   col1  col2
0     1     1
1     1     1
2     1     1
3     1     1
4     1     1
5     1     0
6     0     0
'''
