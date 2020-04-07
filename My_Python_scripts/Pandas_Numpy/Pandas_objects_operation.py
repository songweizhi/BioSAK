import numpy as np
import pandas as pd


# loop over Pandas Series
s = pd.Series(data=np.arange(3), index=['A', 'B', 'C'])
for index, value in s.items():
    print('%s\t%s' % (index, value))

