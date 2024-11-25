import math
import pandas as pd


d = {'col1': [0, 1], 'col2': [10, 100]}

df = pd.DataFrame(data=d)
print(df)


df_log10 = df.map(math.log10)
print(df_log10)

