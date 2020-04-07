# https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
# https://campus.datacamp.com/courses/practicing-statistics-interview-questions-in-python/statistical-experiments-and-significance-testing?ex=13

from statsmodels.stats.multitest import multipletests

p_values = [0.01, 0.05, 0.10, 0.50, 0.99]

p_values_adjusted = multipletests(p_values, alpha=0.1, method='bonferroni')[1]

print(p_values_adjusted)

