__author__ = 'weizhisong'

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.DataFrame(np.random.rand(50, 4), columns=['a', 'b', 'c', 'd'])
ax = df.plot(kind='scatter', x='a', y='b',color='Blue', label='Group 1');
df.plot(kind='scatter', x='c', y='d',color='red', label='Group 2', ax=ax);

plt.show()