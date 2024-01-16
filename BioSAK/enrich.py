import os
import glob
import argparse
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


enrich_usage = '''
=================================== enrich example commands ===================================

BioSAK enrich -i annotation_files -x txt -g grouping.txt -o output_dir -f

# The number of groups needs to be Two!!!

# Example input files:
https://github.com/songweizhi/BioSAK/tree/master/demo_data/enrich

# How it works (https://doi.org/10.1038/s41396-020-00815-8):
Functions that are enriched in the genomes in either group are identified using Mann–Whitney 
U tests followed by a Bonferroni correction with a p value cut-off of 0.05 being considered 
significant. Only significantly different functions with greater than 2-fold mean differences 
are considered to be enriched. Functions detected only in the genomes from one group type are 
considered to be enriched if they existed in at least 50 percent of the genomes in the group.

===============================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def remove_0_from_Pandas_Series(Pandas_Series):
    no_0_num_list = []
    for index, value in Pandas_Series.items():
        if value > 0:
            no_0_num_list.append(value)

    return no_0_num_list


def summarize_stats(output_test, summary_txt):

    sample_1_id = ''
    sample_2_id = ''
    line_num_index = 0
    for ko in open(output_test):
        if line_num_index == 0:
            sample_1_id = ko.strip().split('\t')[1]
            sample_2_id = ko.strip().split('\t')[3]
        line_num_index += 1

    summary_txt_handle = open(summary_txt, 'w')
    summary_txt_handle.write('ID\tP_value\t%s\t%s\tMean_diff\tEnriched_in\n' % (sample_1_id, sample_2_id))
    line_num_index = 0
    for ko in open(output_test):
        if line_num_index > 0:
            ko_split = ko.strip().split('\t')
            ko_id = ko_split[0]
            sample_1_mean = float(ko_split[1])
            sample_1_dectected_pct = float(ko_split[2])
            sample_2_mean = float(ko_split[3])
            sample_2_dectected_pct = float(ko_split[4])
            P_value_adjusted = float(ko_split[6])

            if P_value_adjusted <= 0.05:

                enriched_in = ''
                if (sample_1_mean > 0) and (sample_2_mean == 0):
                    if sample_1_dectected_pct >= 50:
                        enriched_in = sample_1_id
                        mean_diff = 'NA'
                elif (sample_1_mean == 0) and (sample_2_mean > 0):
                    if sample_2_dectected_pct >= 50:
                        enriched_in = sample_2_id
                        mean_diff = 'NA'
                elif (sample_1_mean > 0) and (sample_2_mean > 0):
                    mean_diff = float("{0:.3f}".format(sample_1_mean / sample_2_mean))
                    if mean_diff >= 2:
                        enriched_in = sample_1_id
                    elif mean_diff <= 0.5:
                        enriched_in = sample_2_id

                if (enriched_in == sample_1_id) or (enriched_in == sample_2_id):
                    summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, enriched_in))
        line_num_index += 1
    summary_txt_handle.close()


def enrich(args):

    annotation_file_dir = args['i']
    file_ext            = args['x']
    grouping_file       = args['g']
    op_dir              = args['o']
    force_create_op_dir = args['f']

    annotation_matrix_file  = '%s/data_matrix.txt'  % op_dir
    stats_op_txt            = '%s/stats_result.txt' % op_dir
    summary_txt             = '%s/summary.txt'      % op_dir

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    os.system('mkdir %s' % op_dir)

    # read in grouping file
    grouping_dict = dict()
    group_id_set = set()
    for each_gnm in open(grouping_file):
        each_gnm_split = each_gnm.strip().split('\t')
        grouping_dict[each_gnm_split[0]] = each_gnm_split[1]
        group_id_set.add(each_gnm_split[1])

    group_id_list = sorted([i for i in group_id_set])

    if len(group_id_list) != 2:
        print('The number of groups needs to be TWO, program exited!')
        exit()

    group_1_id = group_id_list[0]
    group_2_id = group_id_list[1]

    annotation_file_re = '%s/*.%s' % (annotation_file_dir, file_ext)
    annotation_file_list = glob.glob(annotation_file_re)
    if len(annotation_file_list) == 0:
        print('Annotation file not found, program exited!')
        exit()

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
            current_MAG_ko_stats_dict_norm[cog] = float("{0:.3f}".format(current_MAG_ko_stats_dict[cog]*100/current_MAG_ko_total))
        ko_dict_of_dict[MAG_id] = current_MAG_ko_stats_dict_norm

    identified_ko_list_sorted = sorted([i for i in identified_ko_set])

    file_out_handle = open(annotation_matrix_file, 'w')
    file_out_handle.write('Group,Genome,%s\n' % (','.join(identified_ko_list_sorted)))
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
    col_index = 0
    for column in column_list:
        if col_index >= 2:
            group1 = df.where(df.Group == group_1_id).dropna()[column]
            group2 = df.where(df.Group == group_2_id).dropna()[column]
            group1_no_zero_pct = len(remove_0_from_Pandas_Series(group1))*100/len(group1)
            group2_no_zero_pct = len(remove_0_from_Pandas_Series(group2))*100/len(group2)
            group1_no_zero_pct = float("{0:.2f}".format(group1_no_zero_pct))
            group2_no_zero_pct = float("{0:.2f}".format(group2_no_zero_pct))
            ko_to_group_detected_pct_dict[column] = [group1_no_zero_pct, group2_no_zero_pct]

            print('Processing %s/%s: %s' % ((col_index - 1), len(column_list)-1, column))

            # store group mean into dict
            group1_mean = float("{0:.2f}".format(sum(group1) / len(group1)))
            group2_mean = float("{0:.2f}".format(sum(group2) / len(group2)))
            ko_to_group_mean_dict[column] = [group1_mean, group2_mean]

            # perform Mann-Whitney U test
            Mann_Whitney_U_test_results[column] = stats.mannwhitneyu(group1, group2)[1]

            ko_id_list.append(column)
            p_value_list.append(stats.ttest_ind(group1, group2, equal_var=False)[1])

        col_index += 1

    p_value_list_adjusted = multipletests(p_value_list, alpha=0.1, method='bonferroni')[1]

    x = 0
    output_test_handle = open(stats_op_txt, 'w')
    output_test_handle.write('KO\t%s\tDectected_pct\t%s\tDectected_pct\tP_value\tP_value_adjusted\n' % (group_1_id, group_2_id))
    while x < len(ko_id_list):
        current_p = float("{0:.3f}".format(p_value_list[x]))
        current_p_adjusted = float("{0:.3f}".format(p_value_list_adjusted[x]))
        current_p_group_mean = [str(i) for i in ko_to_group_mean_dict[ko_id_list[x]]]
        current_p_detected_pct = [str(i) for i in ko_to_group_detected_pct_dict[ko_id_list[x]]]

        group_1_mean = current_p_group_mean[0]
        group_2_mean = current_p_group_mean[1]
        group_1_no_zero_pct = current_p_detected_pct[0]
        group_2_no_zero_pct = current_p_detected_pct[1]

        output_test_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id_list[x], group_1_mean, group_1_no_zero_pct, group_2_mean, group_2_no_zero_pct, current_p, current_p_adjusted))
        x += 1
    output_test_handle.close()

    # summarize stats
    summarize_stats(stats_op_txt, summary_txt)

    # file report
    print('Results of enrichment analysis exported to: %s' % summary_txt)
    print('Done!')


if __name__ == '__main__':

    enrich_parser = argparse.ArgumentParser(usage=enrich_usage)
    enrich_parser.add_argument('-i', required=True,                         help='annotation files')
    enrich_parser.add_argument('-x', required=True,                         help='file extension')
    enrich_parser.add_argument('-g', required=True,                         help='grouping file')
    enrich_parser.add_argument('-o', required=True,                         help='output directory')
    enrich_parser.add_argument('-f', required=False, action="store_true",   help='force overwrite')
    args = vars(enrich_parser.parse_args())
    enrich(args)
