import os
import glob
import argparse
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


enrich_usage = '''
======================================= enrich example commands =======================================

BioSAK enrich -i annotation_files -x txt -g grouping.txt -o output_dir
BioSAK enrich -i annotation_files -x txt -g grouping.txt -o output_dir -bc

# Note:
1. The number of genome groups has to be TWO!!!
2. Group name should NOT be 1 and 0. use some alphabet or words instead.

# Example input files: https://github.com/songweizhi/BioSAK/tree/master/demo_data/enrich

# please refers to https://doi.org/10.1038/s41396-020-00815-8 for how it works:
Functions that are enriched in the genomes in either group are identified using Mann–Whitney U 
tests followed by a Bonferroni correction with a p value cut-off of 0.05 being considered 
significant. Only significantly different functions with greater than, by default, 2-fold mean 
differences are considered to be enriched. Functions detected only in the genomes from one group 
type are considered to be enriched if they existed in at least 50 percent of the genomes in the group.

=======================================================================================================
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


def summarize_stats(output_test, fold_diff_cutoff, ko_desc_dict, fun_to_note_dict, file_prefix, output_dir):

    fold_diff_big = fold_diff_cutoff
    fold_diff_small = 1/fold_diff_cutoff
    if (1/fold_diff_cutoff) > fold_diff_cutoff:
        fold_diff_big = 1/fold_diff_cutoff
        fold_diff_small = fold_diff_cutoff

    sample_1_id = ''
    sample_2_id = ''
    line_num_index = 0
    for ko in open(output_test):
        if line_num_index == 0:
            sample_1_id = ko.strip().split('\t')[1]
            sample_2_id = ko.strip().split('\t')[3]
        line_num_index += 1

    summary_txt_sample_1 = '%s/enriched_in_%s.txt' % (output_dir, sample_1_id)
    summary_txt_sample_2 = '%s/enriched_in_%s.txt' % (output_dir, sample_2_id)
    if file_prefix != '':
        summary_txt_sample_1 = '%s/%s_enriched_in_%s.txt' % (output_dir, file_prefix, sample_1_id)
        summary_txt_sample_2 = '%s/%s_enriched_in_%s.txt' % (output_dir, file_prefix, sample_2_id)

    summary_txt_sample_1_handle = open(summary_txt_sample_1, 'w')
    summary_txt_sample_2_handle = open(summary_txt_sample_2, 'w')
    if len(fun_to_note_dict) == 0:
        summary_txt_sample_1_handle.write('ID\tP_value\t%s\t%s\tMean_diff\tDescription\n' % (sample_1_id, sample_2_id))
        summary_txt_sample_2_handle.write('ID\tP_value\t%s\t%s\tMean_diff\tDescription\n' % (sample_1_id, sample_2_id))
    else:
        summary_txt_sample_1_handle.write('ID\tP_value\t%s\t%s\tMean_diff\tDescription\tNote\n' % (sample_1_id, sample_2_id))
        summary_txt_sample_2_handle.write('ID\tP_value\t%s\t%s\tMean_diff\tDescription\tNote\n' % (sample_1_id, sample_2_id))
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

                mean_diff = 'NA'
                enriched_in = ''
                if (sample_1_mean > 0) and (sample_2_mean == 0):
                    if sample_1_dectected_pct >= 50:
                        enriched_in = sample_1_id
                elif (sample_1_mean == 0) and (sample_2_mean > 0):
                    if sample_2_dectected_pct >= 50:
                        enriched_in = sample_2_id
                elif (sample_1_mean > 0) and (sample_2_mean > 0):
                    mean_diff = float("{0:.3f}".format(sample_1_mean/sample_2_mean))
                    if mean_diff >= fold_diff_big:
                        enriched_in = sample_1_id
                    elif mean_diff <= fold_diff_small:
                        enriched_in = sample_2_id

                if mean_diff != 'NA':
                    if mean_diff < 1:
                        mean_diff = 1/mean_diff
                        mean_diff = float("{0:.3f}".format(mean_diff))

                if enriched_in == sample_1_id:
                    if len(fun_to_note_dict) == 0:
                        summary_txt_sample_1_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, ko_desc_dict.get(ko_id, 'na')))
                    else:
                        summary_txt_sample_1_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, ko_desc_dict.get(ko_id, 'na'), fun_to_note_dict.get(ko_id, 'na')))
                if enriched_in == sample_2_id:
                    if len(fun_to_note_dict) == 0:
                        summary_txt_sample_2_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, ko_desc_dict.get(ko_id, 'na')))
                    else:
                        summary_txt_sample_2_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id, P_value_adjusted, sample_1_mean, sample_2_mean, mean_diff, ko_desc_dict.get(ko_id, 'na'), fun_to_note_dict.get(ko_id, 'na')))
        line_num_index += 1
    summary_txt_sample_1_handle.close()
    summary_txt_sample_2_handle.close()


def enrich(args):

    op_prefix           = args['p']
    annotation_file_dir = args['i']
    file_ext            = args['x']
    grouping_file       = args['g']
    op_dir              = args['o']
    fold_diff_cutoff    = args['diff']
    perform_bc          = args['bc']
    force_create_op_dir = args['f']

    ####################################################################################################################

    # read in grouping file
    group_id_set = set()
    grouping_dict = dict()
    gnm_with_grouping_set = set()
    for each_gnm in open(grouping_file):
        each_gnm_split = each_gnm.strip().split('\t')
        grouping_dict[each_gnm_split[0]] = each_gnm_split[1]
        group_id_set.add(each_gnm_split[1])
        gnm_with_grouping_set.add(each_gnm_split[0])

    if len(group_id_set) != 2:
        print('The number of groups needs to be TWO, program exited!')
        exit()

    # get annotation file list
    annotation_file_re = '%s/*.%s' % (annotation_file_dir, file_ext)
    annotation_file_list = glob.glob(annotation_file_re)
    annotation_file_list_id = [os.path.splitext(os.path.basename(i))[0] for i in annotation_file_list]
    if len(annotation_file_list_id) == 0:
        print('Annotation file not found, program exited!')
        exit()

    # get shared genomes
    shared_gnm_set = set(gnm_with_grouping_set).intersection(annotation_file_list_id)
    if len(shared_gnm_set) == 0:
        print('No genome shared by %s and %s, program exited!' % (grouping_file, annotation_file_dir))
        exit()

    gnms_uniq_to_grouping_file = set()
    for g1 in gnm_with_grouping_set:
        if g1 not in shared_gnm_set:
            gnms_uniq_to_grouping_file.add(g1)

    gnms_uniq_to_annotation_dir = set()
    for g2 in annotation_file_list_id:
        if g2 not in shared_gnm_set:
            gnms_uniq_to_annotation_dir.add(g2)

    # report
    if len(gnms_uniq_to_grouping_file) > 0:
        print('Genomes uniq to %s (%s): %s' % (grouping_file, len(gnms_uniq_to_grouping_file), ','.join(sorted(list(gnms_uniq_to_grouping_file)))))
    if len(gnms_uniq_to_annotation_dir) > 0:
        print('Genomes uniq to %s (%s): %s' % (annotation_file_dir, len(gnms_uniq_to_annotation_dir), ','.join(sorted(list(gnms_uniq_to_annotation_dir)))))

    print('%s genomes found in both %s and %s.' % (len(shared_gnm_set), grouping_file, annotation_file_dir))
    print('Will perform enrichment analysis on the %s genomes.' % len(shared_gnm_set))

    ####################################################################################################################

    annotation_matrix_file  = '%s/data_matrix.txt'  % op_dir
    stats_op_txt            = '%s/stats_result.txt' % op_dir
    if op_prefix != '':
        annotation_matrix_file  = '%s/%s_data_matrix.txt'  % (op_dir, op_prefix)
        stats_op_txt            = '%s/%s_stats_result.txt' % (op_dir, op_prefix)

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    ####################################################################################################################

    annotation_file_list_only_shared = set()
    for each_file in annotation_file_list:
        _, f_base, _ = sep_path_basename_ext(each_file)
        if f_base in shared_gnm_set:
            annotation_file_list_only_shared.add(each_file)

    group_id_list = sorted([i for i in group_id_set])
    group_1_id = group_id_list[0]
    group_2_id = group_id_list[1]

    identified_ko_set = set()
    ko_dict_of_dict = {}
    ko_desc_dict = dict()
    fun_id_to_note_dict = dict()
    for annotation_file in annotation_file_list_only_shared:
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
                if len(each_line_split) >= 3:
                    ko_desc_dict[each_line_split[0]] = each_line_split[2]
                if len(each_line_split) == 4:
                    fun_id_to_note_dict[each_line_split[0]] = each_line_split[3]
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
        if col_index >= 1:
            print('Processing %s/%s: %s' % ((col_index - 1), len(column_list)-1, column))
            group1 = df.where(df.Group == group_1_id).dropna()[column]
            group2 = df.where(df.Group == group_2_id).dropna()[column]

            # to calculate the presence coverage, nothing about the enrichment stats
            group1_no_zero_pct = len(remove_0_from_Pandas_Series(group1))*100/len(group1)
            group2_no_zero_pct = len(remove_0_from_Pandas_Series(group2))*100/len(group2)
            group1_no_zero_pct = float("{0:.2f}".format(group1_no_zero_pct))
            group2_no_zero_pct = float("{0:.2f}".format(group2_no_zero_pct))
            ko_to_group_detected_pct_dict[column] = [group1_no_zero_pct, group2_no_zero_pct]

            # store group mean into dict
            group1_mean = float("{0:.2f}".format(sum(group1) / len(group1)))
            group2_mean = float("{0:.2f}".format(sum(group2) / len(group2)))
            ko_to_group_mean_dict[column] = [group1_mean, group2_mean]

            # perform Mann-Whitney U test
            Mann_Whitney_U_test_results[column] = stats.mannwhitneyu(group1, group2)[1]

            ko_id_list.append(column)
            p_value_list.append(stats.ttest_ind(group1, group2, equal_var=False)[1])

        col_index += 1

    if perform_bc is True:
        p_value_list_adjusted = multipletests(p_value_list, alpha=0.1, method='bonferroni')[1]
    else:
        p_value_list_adjusted = p_value_list

    x = 0
    output_test_handle = open(stats_op_txt, 'w')
    output_test_handle.write('KO\t%s\tDectected_pct\t%s\tDectected_pct\tP_value\tP_value_adjusted\tDescription\n' % (group_1_id, group_2_id))
    while x < len(ko_id_list):
        current_p = float("{0:.3f}".format(p_value_list[x]))
        current_p_adjusted = float("{0:.3f}".format(p_value_list_adjusted[x]))
        current_p_group_mean = [str(i) for i in ko_to_group_mean_dict[ko_id_list[x]]]
        current_p_detected_pct = [str(i) for i in ko_to_group_detected_pct_dict[ko_id_list[x]]]
        group_1_mean = current_p_group_mean[0]
        group_2_mean = current_p_group_mean[1]
        group_1_no_zero_pct = current_p_detected_pct[0]
        group_2_no_zero_pct = current_p_detected_pct[1]
        output_test_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ko_id_list[x], group_1_mean, group_1_no_zero_pct, group_2_mean, group_2_no_zero_pct, current_p, current_p_adjusted, ko_desc_dict.get(ko_id_list[x], 'na')))
        x += 1
    output_test_handle.close()

    # summarize stats
    summarize_stats(stats_op_txt, fold_diff_cutoff, ko_desc_dict, fun_id_to_note_dict, op_prefix, op_dir)

    # file report
    print('Done!')


if __name__ == '__main__':

    enrich_parser = argparse.ArgumentParser(usage=enrich_usage)
    enrich_parser.add_argument('-p',    required=False, default='',             help='prefix of output files')
    enrich_parser.add_argument('-i',    required=True,                          help='annotation files')
    enrich_parser.add_argument('-x',    required=True,                          help='file extension')
    enrich_parser.add_argument('-g',    required=True,                          help='grouping file')
    enrich_parser.add_argument('-o',    required=True,                          help='output directory')
    enrich_parser.add_argument('-diff', required=False, default=2, type=float,  help='minimum fold difference, default is 2')
    enrich_parser.add_argument('-bc',   required=False, action="store_true",    help='perform Bonferroni correction')
    enrich_parser.add_argument('-f',    required=False, action="store_true",    help='force overwrite')
    args = vars(enrich_parser.parse_args())
    enrich(args)
