import scipy
from scipy import stats
from matplotlib import pyplot


group1 = [1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,8,8,9,0, 10]
group1 = [1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,7,7,7,8,8,9,0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]


sample_normality = 'No'
if float(stats.normaltest(group1).pvalue) > 0.05:
    sample_normality = 'Yes'

pyplot.hist(group1)
pyplot.show()
print('Sample normality: %s' % sample_normality)


print('\nstats.normaltest')
print(stats.normaltest(group1))

print('\nscipy.stats.shapiro')
print(scipy.stats.shapiro(group1))

