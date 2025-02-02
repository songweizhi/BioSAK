import os
import glob
import math
import argparse
import pandas as pd
from functools import reduce


get_abd3_stats_usage = '''
===================== get_abd3_stats example commands =====================

BioSAK get_abd3_stats -i get_abd2_mapping_op_dir -o op_dir -f 
BioSAK get_abd3_stats -i cov_rpkm_stat_files -o get_abd3_stats_op_dir -f 

# Note
0 becomes na after taking log.

===========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_abd3_stats(args):

    rpkm_stat_cov_dir   = args['i']
    op_dir              = args['o']
    force_overwrite     = args['f']

    # check input files
    if os.path.isdir(rpkm_stat_cov_dir) is False:
        print('%s not found, program exited!' % rpkm_stat_cov_dir)
        exit()

    rpkm_file_re    = '%s/*.rpkm' % rpkm_stat_cov_dir
    stat_file_re    = '%s/*.stat' % rpkm_stat_cov_dir
    rpkm_file_list  = glob.glob(rpkm_file_re)
    stat_file_list  = glob.glob(stat_file_re)

    if len(rpkm_file_list) == 0:
        print('.rpkm file not found, program exited!')
        exit()

    if len(stat_file_list) == 0:
        print('.stat file not found, program exited!')
        exit()

    # define output file name
    gnm_level_rpkm_dir       = '%s/genome_level_rpkm'   % op_dir
    combined_rpkm_file       = '%s/combined.rpkm'       % op_dir
    combined_rpkm_file_log10 = '%s/combined.log10.rpkm' % op_dir

    # create op dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output directory detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gnm_level_rpkm_dir)

    # read in stat file
    sample_total_read_num_dict = dict()
    for stat_file in sorted(stat_file_list):
        _, _, stat_base, _ = sep_path_basename_ext(stat_file)
        stat_file_lines = open(stat_file).readlines()
        info_line_split = stat_file_lines[-1].strip().split(' ')
        while '' in info_line_split:
            info_line_split.remove('')
        read_num = int(info_line_split[3].replace(',', ''))
        sample_total_read_num_dict[stat_base] = read_num

    # summarise rpkm on genome level
    sample_id_list = []
    for rpkm_file in sorted(rpkm_file_list):
        rpkm_name, rpkm_path, rpkm_base, rpkm_ext = sep_path_basename_ext(rpkm_file)
        sample_id = rpkm_base.split('.sorted_filtered')[0]
        sample_id_list.append(sample_id)
        op_txt = '%s/%s.cal_bin.rpkm' % (gnm_level_rpkm_dir, sample_id)
        rpkm_df = pd.read_csv(rpkm_file, sep='\t', header=4)
        subset_df_col_list = ['#Name', 'Length', 'Reads']
        sample_total_read_num = sample_total_read_num_dict[sample_id]
        subset_df = rpkm_df.loc[:, subset_df_col_list]
        subset_df['bin_name'] = subset_df['#Name'].str.rsplit('_', n=1, expand=True)[0]
        subset_df['Length'] = subset_df['Length'].astype(int)
        subset_df['Reads'] = subset_df['Reads'].astype(int)

        # calulate RPKM
        op_df = subset_df.groupby(["bin_name"])[['Length', 'Reads']].sum()
        op_df[sample_id] = (op_df['Reads']*1000000000)/(op_df['Length']*sample_total_read_num*2)
        op_df.to_csv(op_txt, sep='\t')

    # combine the results
    rpkm_df_list = []
    for sample_id in sample_id_list:
        gnm_level_rpkm_file = '%s/%s.cal_bin.rpkm' % (gnm_level_rpkm_dir, sample_id)
        gnm_level_rpkm_df = pd.read_csv(gnm_level_rpkm_file, sep='\t')
        gnm_level_rpkm_df = gnm_level_rpkm_df.drop(["Length", "Reads"], axis=1)
        rpkm_df_list.append(gnm_level_rpkm_df)

    # write out combined_rpkm_file
    rpkm_df_combined = reduce(lambda left, right: pd.merge(left, right, on='bin_name',how='outer'), rpkm_df_list)
    rpkm_df_combined.to_csv(combined_rpkm_file, sep='\t',index=False)

    # write out combined_rpkm_file_log10
    combined_rpkm_file_log10_handle = open(combined_rpkm_file_log10, 'w')
    for each_line in open(combined_rpkm_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('bin_name'):
            combined_rpkm_file_log10_handle.write(each_line)
        else:
            value_list = [each_line_split[0]]
            for each_value in each_line_split[1:]:
                each_value = float(each_value)
                if each_value == 0:
                    value_list.append('na')
                else:
                    each_value_log = math.log(each_value)
                    value_list.append(each_value_log)
            combined_rpkm_file_log10_handle.write('%s\n' % '\t'.join([str(i) for i in value_list]))
    combined_rpkm_file_log10_handle.close()

    # final report
    print('Done!')


if __name__ == '__main__':

    get_abd3_stats_parser = argparse.ArgumentParser()
    get_abd3_stats_parser.add_argument('-i', required=True,                          help='output directory from get_abd2_mapping')
    get_abd3_stats_parser.add_argument('-o', required=True,                          help='output directory')
    get_abd3_stats_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    args = vars(get_abd3_stats_parser.parse_args())
    get_abd3_stats(args)
