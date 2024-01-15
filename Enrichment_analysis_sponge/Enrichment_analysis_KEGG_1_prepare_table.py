import os
import glob
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def remove_0_from_Pandas_Series(Pandas_Series):

    no_0_num_list = []
    for index, value in Pandas_Series.items():
        if value > 0:
            no_0_num_list.append(value)

    return no_0_num_list


def enrich(args):

    ########################################################################################################################

    # file in
    annotation_file_dir     = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber'
    file_ext                = 'txt'
    grouping_file           = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/grouping.txt'

    # file out
    annotation_matrix_file  = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber.csv'
    stats_op_txt            = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber_Mann_Whitney_U.tab'

    ########################################################################################################################

    # read in grouping file
    grouping_dict = dict()
    for each_gnm in open(grouping_file):
        each_gnm_split = each_gnm.strip().split('\t')
        grouping_dict[each_gnm_split[0]] = each_gnm_split[1]

    annotation_file_re = '%s/*.%s' % (annotation_file_dir, file_ext)
    annotation_file_list = glob.glob(annotation_file_re)

    identified_ko_set = set()
    ko_dict_of_dict = {}
    for annotation_file in annotation_file_list:
        f_path, MAG_id, f_ext = sep_path_basename_ext(annotation_file)
        current_MAG_ko_stats_dict = {}
        current_MAG_ko_total = 0
        line_num_index = 0
        for each_line in open(annotation_file):
            if line_num_index > 0:
                each_line_split = each_line.strip().split('\t')
                identified_ko_set.add(each_line_split[0])
                current_MAG_ko_stats_dict[each_line_split[0]] = int(each_line_split[1])
                current_MAG_ko_total += int(each_line_split[1])
            line_num_index += 1

        current_MAG_ko_stats_dict_norm = {}
        for cog in current_MAG_ko_stats_dict:
            current_MAG_ko_stats_dict_norm[cog] = float("{0:.3f}".format(current_MAG_ko_stats_dict[cog] * 100 / current_MAG_ko_total))
        ko_dict_of_dict[MAG_id] = current_MAG_ko_stats_dict_norm

    identified_ko_list_sorted = sorted([i for i in identified_ko_set])

    file_out_handle = open(annotation_matrix_file, 'w')
    file_out_handle.write('Source,MAG,%s\n' % (','.join(identified_ko_list_sorted)))
    for MAG in ko_dict_of_dict:
        current_MAG_ko_stats_list = []
        current_MAG_ko_stats = ko_dict_of_dict[MAG]
        for ko_id in identified_ko_list_sorted:
            ko_id_num = 0
            if ko_id in current_MAG_ko_stats:
                ko_id_num = current_MAG_ko_stats[ko_id]
            current_MAG_ko_stats_list.append(ko_id_num)

        current_MAG_ko_stats_list_str = [str(i) for i in current_MAG_ko_stats_list]

        if MAG not in grouping_dict:
            print('ID %s not found in %s, program exited!' % (MAG, grouping_file))
            print('Program exited!')
            exit()

        mag_group = grouping_dict[MAG]
        file_out_handle.write('%s,%s,%s\n' % (mag_group, MAG, ','.join(current_MAG_ko_stats_list_str)))
    file_out_handle.close()

    ########################################################################################################################

    df = pd.read_csv(annotation_matrix_file)

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
    output_test_handle = open(stats_op_txt, 'w')
    output_test_handle.write('KO\tsponge\tDectected_pct\tseawater\tDectected_pct\tP_value\tP_value_adjusted\n')
    while x < len(ko_id_list):

        current_p               = float("{0:.3f}".format(p_value_list[x]))
        current_p_adjusted      = float("{0:.3f}".format(p_value_list_adjusted[x]))
        current_p_group_mean    = [str(i) for i in ko_to_group_mean_dict[ko_id_list[x]]]
        current_p_detected_pct  = [str(i) for i in ko_to_group_detected_pct_dict[ko_id_list[x]]]

        output_test_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id_list[x], current_p_group_mean[0], current_p_detected_pct[0], current_p_group_mean[1], current_p_detected_pct[1], current_p, current_p_adjusted))

        x += 1

    output_test_handle.close()

