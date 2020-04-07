from scipy import stats


# normality test
print(stats.normaltest([1,2,3,4,5,6,7,8,9,6,4,5,6,7,4,3,5,6,7,8,6,4,5,7,8,9,4,2,4,6,7,4,6,6]))
# NormaltestResult(statistic=0.4034400613897061, pvalue=0.8173237214675394)

print(stats.normaltest([1,1,1,1,1,11,1,1,1,6,6,6,6,6,69,9,9,9,9,9,9,9,9,9,9]))
# NormaltestResult(statistic=55.677348367888854, pvalue=8.124888661957071e-13)

