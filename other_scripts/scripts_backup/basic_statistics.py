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

