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
csv_file =      '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber.csv'
output_test =   '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber_Mann_Whitney_U.tab'


# df = pd.read_csv(csv_file, index_col=0)
df = pd.read_csv(csv_file)

# get a list of all columns in the dataframe without the Group column
column_list = [x for x in df.columns if x != 'Group']

# loop over column_list and execute code explained above
Mann_Whitney_U_test_results = {}
ko_to_group_mean_dict = {}
ko_to_group_detected_pct_dict = {}
p_value_list = []
ko_id_list = []
n = 0
for column in column_list:
    if n >= 2:
        group1 = df.where(df.Source == 'sponge').dropna()[column]
        group2 = df.where(df.Source == 'seawater').dropna()[column]
        group1_no_zero_pct = len(remove_0_from_Pandas_Series(group1)) * 100 / len(group1)
        group2_no_zero_pct = len(remove_0_from_Pandas_Series(group2)) * 100 / len(group2)
        group1_no_zero_pct = float("{0:.2f}".format(group1_no_zero_pct))
        group2_no_zero_pct = float("{0:.2f}".format(group2_no_zero_pct))
        ko_to_group_detected_pct_dict[column] = [group1_no_zero_pct, group2_no_zero_pct]

        print('Processing %s/%s: %s' % (n, len(column_list), column))

        # store group mean into dict
        group1_mean = float("{0:.2f}".format(sum(group1) / len(group1)))
        group2_mean = float("{0:.2f}".format(sum(group2) / len(group2)))
        ko_to_group_mean_dict[column] = [group1_mean, group2_mean]

        # perform Mann-Whitney U test
        Mann_Whitney_U_test_results[column] = stats.mannwhitneyu(group1, group2)[1]

        ko_id_list.append(column)
        p_value_list.append(stats.ttest_ind(group1, group2, equal_var=False)[1])

    n += 1

p_value_list_adjusted = multipletests(p_value_list, alpha=0.1, method='bonferroni')[1]


x = 0
output_test_handle = open(output_test, 'w')
output_test_handle.write('KO\tsponge\tDectected_pct\tseawater\tDectected_pct\tP_value\tP_value_adjusted\n')
while x < len(ko_id_list):

    current_p               = float("{0:.3f}".format(p_value_list[x]))
    current_p_adjusted      = float("{0:.3f}".format(p_value_list_adjusted[x]))
    current_p_group_mean    = [str(i) for i in ko_to_group_mean_dict[ko_id_list[x]]]
    current_p_detected_pct  = [str(i) for i in ko_to_group_detected_pct_dict[ko_id_list[x]]]

    output_test_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id_list[x], current_p_group_mean[0], current_p_detected_pct[0], current_p_group_mean[1], current_p_detected_pct[1], current_p, current_p_adjusted))

    x += 1

output_test_handle.close()

