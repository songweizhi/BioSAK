import scipy
from scipy import stats
import pandas as pd
from pandas import DataFrame
from matplotlib import pyplot
from statsmodels.stats.multitest import multipletests


def remove_0_from_Pandas_Series(Pandas_Series):

    no_0_num_list = []
    for index, value in Pandas_Series.items():
        if value > 0:
            no_0_num_list.append(value)

    return no_0_num_list


# file in
#csv_file =     '/Users/songweizhi/Desktop/COG_enrichment_analysis/Herbivore_specialisation.csv'
csv_file =      '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_analysis.csv'
output_test =   '/Users/songweizhi/Desktop/COG_enrichment_analysis/COG_enrichment_test_results.tab'


# df = pd.read_csv(csv_file, index_col=0)
df = pd.read_csv(csv_file)

# get a list of all columns in the dataframe without the Group column
column_list = [x for x in df.columns if x != 'Group']

# loop over column_list and execute code explained above
t_test_results = {}
cog_to_group_mean_dict = {}
cog_to_group_mean_dict_no_zero = {}
n = 0
a = 0
b = 0
p_value_list = []
COG_id_list = []
for column in column_list:

    if n >= 2:

        group1 = df.where(df.Source == 'kelp-associated').dropna()[column]
        group2 = df.where(df.Source == 'planktonic').dropna()[column]

        group1_no_zero = remove_0_from_Pandas_Series(group1)
        group2_no_zero = remove_0_from_Pandas_Series(group2)

        group1_normality = 'No'
        if len(group1) >= 20:
            if float(stats.normaltest(group1).pvalue) > 0.05:
                group1_normality = 'Yes'

        group2_normality = 'No'
        if len(group2) >= 20:
            if float(stats.normaltest(group2).pvalue) > 0.05:
                group2_normality = 'Yes'

        group1_no_zero_normality = 'No'
        if len(group1_no_zero) >= 20:
            if float(stats.normaltest(group1_no_zero).pvalue) > 0.05:
                group1_no_zero_normality = 'Yes'

        group2_no_zero_normality = 'No'
        if len(group2_no_zero) >= 20:
            if float(stats.normaltest(group2_no_zero).pvalue) > 0.05:
                group2_no_zero_normality = 'Yes'

        if (group1_normality == 'Yes') and (group2_normality == 'Yes'):
            a += 1

        if (group1_no_zero_normality == 'Yes') and (group2_no_zero_normality == 'Yes'):

            print('Processing the %sth column: %s' % (n, column))

            group1_mean = float("{0:.2f}".format(sum(group1) / len(group1)))
            group2_mean = float("{0:.2f}".format(sum(group2) / len(group2)))

            group1_no_zero_mean = float("{0:.2f}".format(sum(group1_no_zero) / len(group1_no_zero)))
            group2_no_zero_mean = float("{0:.2f}".format(sum(group2_no_zero) / len(group2_no_zero)))

            cog_to_group_mean_dict[column] = [group1_mean, group2_mean]
            cog_to_group_mean_dict_no_zero[column] = [group1_no_zero_mean, group2_no_zero_mean]

            t_test_results[column] = stats.ttest_ind(group1, group2, equal_var=False)[1]  # If equal_var=False, perform Welch’s t-test, which does not assume equal population variance [2].

            COG_id_list.append(column)
            p_value_list.append(stats.ttest_ind(group1, group2, equal_var=False)[1])

            b += 1

        # print('Normality Test: %s' % group1_normality)
        # pyplot.hist(group1)
        # pyplot.show()
        #
        # print('Normality Test: %s' % group1_no_zero_normality)
        # pyplot.hist(group1_no_zero)
        # pyplot.show()

    n += 1

p_value_list_adjusted = multipletests(p_value_list, alpha=0.1, method='bonferroni')[1]

print('a')
print(a)
print('b')
print(b)

x = 0
output_test_handle = open(output_test, 'w')
output_test_handle.write('COG\tKelp\tPlanktonic\tP_value\tP_value_adjusted\n')
while x < len(COG_id_list):

    current_p            = float("{0:.3f}".format(p_value_list[x]))
    current_p_adjusted   = float("{0:.3f}".format(p_value_list_adjusted[x]))
    current_p_group_mean = [str(i) for i in cog_to_group_mean_dict_no_zero[COG_id_list[x]]]

    if current_p_adjusted <= 0.05:
        output_test_handle.write('%s\t%s\t%s\t%s\n' % (COG_id_list[x], '\t'.join(current_p_group_mean), current_p, current_p_adjusted))

    x += 1
output_test_handle.close()


'''
overall:


normality:
0

no 0 normality:
43

'''
