import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


'''

Python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_KEGG.py -in Kelp_dRep_fun_stats -out Kelp_dRep_fun_stats.txt -in_percent
Python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_KEGG.py -in Tara_NM_fun_stats -out Tara_NM_fun_stats.txt -in_percent

Rscript ~/PycharmProjects/MyBioTools/MyBioTools/COG_boxplot_last1row.R -i Kelp_dRep_fun_stats.txt -o Kelp_dRep_fun_stats.png
Rscript ~/PycharmProjects/MyBioTools/MyBioTools/COG_boxplot_last1row.R -i ko_B_percent_df_no_NA.csv -o ko_B_percent_df_no_NA.png

'''


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


def get_bin_to_ko_percent_dict(pwd_ko_stats, ko_id_all):

    bin_to_ko_percent = {}
    for each_ko in open(pwd_ko_stats):
        if not each_ko.startswith('Level'):
            each_ko_split = each_ko.strip().split('\t')
            each_ko_id = each_ko_split[1]
            each_ko_percent = float(each_ko_split[3])
            if each_ko_id not in ko_id_all:
                ko_id_all.append(each_ko_id)

            bin_to_ko_percent[each_ko_id] = each_ko_percent

    return bin_to_ko_percent


# input files
KEGG_annot_results_folder = '/Users/songweizhi/Desktop/SpongeMAGs_dRep_KEGG_wd'

# output files
ko_B_percent_csv_no_NA = '/Users/songweizhi/Desktop/ko_B_percent_df_no_NA.csv'
ko_C_percent_csv_no_NA = '/Users/songweizhi/Desktop/ko_C_percent_df_no_NA.csv'


# DB file
KO_description_B = '/Users/songweizhi/DB/KEGG_DB/KO_description_B.txt'
KO_description_C = '/Users/songweizhi/DB/KEGG_DB/KO_description_C.txt'


columns_to_drop_B = []
columns_to_drop_C = []


KO_description_B_dict = {}
for each_ko_B in open(KO_description_B):
    each_ko_B_split = each_ko_B.strip().split('\t')
    KO_description_B_dict[each_ko_B_split[0]] = each_ko_B_split[1]


KO_description_C_dict = {}
for each_ko_C in open(KO_description_C):
    each_ko_C_split = each_ko_C.strip().split('\t')
    KO_description_C_dict[each_ko_C_split[0]] = each_ko_C_split[1]


KEGG_wd_folder_list = get_no_hidden_folder_list(KEGG_annot_results_folder)

bin_id_list = []
ko_id_B_all = []
ko_id_C_all = []
ko_percent_all_bins_B = {}
ko_percent_all_bins_C = {}
for each_folder in KEGG_wd_folder_list:

    bin_id = each_folder.split('_KEGG_wd')[0]
    bin_id_list.append(bin_id)

    pwd_ko_stats_A = '%s/%s/%s_ko_stats_A.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_B = '%s/%s/%s_ko_stats_B.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_C = '%s/%s/%s_ko_stats_C.txt' % (KEGG_annot_results_folder, each_folder, bin_id)
    pwd_ko_stats_D = '%s/%s/%s_ko_stats_D.txt' % (KEGG_annot_results_folder, each_folder, bin_id)

    bin_to_ko_percent_B = get_bin_to_ko_percent_dict(pwd_ko_stats_B, ko_id_B_all)
    bin_to_ko_percent_C = get_bin_to_ko_percent_dict(pwd_ko_stats_C, ko_id_C_all)

    ko_percent_all_bins_B[bin_id] = bin_to_ko_percent_B
    ko_percent_all_bins_C[bin_id] = bin_to_ko_percent_C


bin_id_list = sorted(bin_id_list)
ko_id_B_all = sorted(ko_id_B_all)
ko_id_C_all = sorted(ko_id_C_all)


KO_description_B_dict['NA'] = 'NA'
KO_description_C_dict['NA'] = 'NA'

ko_id_B_all_description = ['B_%s_%s' % (i, KO_description_B_dict[i]) for i in ko_id_B_all]
ko_id_C_all_description = ['C_%s_%s' % (i, KO_description_C_dict[i]) for i in ko_id_C_all]


ko_percent_lol_B = []
ko_percent_lol_C = []
for each_bin in bin_id_list:
    current_bin_ko_B_percent_dict = ko_percent_all_bins_B[each_bin]
    current_bin_ko_C_percent_dict = ko_percent_all_bins_C[each_bin]
    percent_list_B = []
    for ko_id_B in ko_id_B_all:
        ko_B_percent = 0
        if ko_id_B in current_bin_ko_B_percent_dict:
            ko_B_percent = current_bin_ko_B_percent_dict[ko_id_B]
        percent_list_B.append(ko_B_percent)

    percent_list_C = []
    for ko_id_C in ko_id_C_all:
        ko_C_percent = 0
        if ko_id_C in current_bin_ko_C_percent_dict:
            ko_C_percent = current_bin_ko_C_percent_dict[ko_id_C]
        percent_list_C.append(ko_C_percent)

    ko_percent_lol_B.append(percent_list_B)
    ko_percent_lol_C.append(percent_list_C)

    # percent_list_str_B = [str(i) for i in percent_list_B]
    # percent_list_str_C = [str(i) for i in percent_list_C]
    # for_out_B = '%s,%s' % (each_bin, ','.join(percent_list_str_B))
    # for_out_C = '%s,%s' % (each_bin, ','.join(percent_list_str_C))
    # print(for_out_B)
    # print(for_out_C)


ko_B_percent_df = pd.DataFrame(np.array(ko_percent_lol_B), index=bin_id_list, columns=ko_id_B_all_description)
ko_C_percent_df = pd.DataFrame(np.array(ko_percent_lol_C), index=bin_id_list, columns=ko_id_C_all_description)


# ko_B_percent_df = pd.DataFrame(np.array(ko_percent_lol_B), index=bin_id_list, columns=ko_id_B_all)
# ko_C_percent_df = pd.DataFrame(np.array(ko_percent_lol_C), index=bin_id_list, columns=ko_id_C_all)


ko_B_percent_df_no_NA = ko_B_percent_df.loc[:, ko_B_percent_df.columns != 'B_NA_NA']
ko_C_percent_df_no_NA = ko_C_percent_df.loc[:, ko_C_percent_df.columns != 'C_NA_NA']


ko_B_percent_df_no_NA.to_csv(ko_B_percent_csv_no_NA, header=True)
ko_C_percent_df_no_NA.to_csv(ko_C_percent_csv_no_NA, header=True)


# print(ko_B_percent_df.shape)
# print(ko_B_percent_df.index)
# for each in ko_B_percent_df.index:
#     print(each)
#print(ko_B_percent_df)
# ko_B_percent_df_NA = ko_B_percent_df.loc[:, ko_B_percent_df.columns == 'NA']
# print(ko_B_percent_df_NA)
# print(ko_percent_df.index)
# print(ko_percent_df.columns)
# print(ko_percent_df.values)
# print(ko_C_percent_df.shape)
# print(ko_C_percent_df)
# ko_C_percent_df_no_NA.boxplot(grid=False, figsize=(100, 30))
# plt.savefig('/Users/songweizhi/Desktop/test_C.png', dpi=300)
# ko_B_percent_df_no_NA.boxplot(grid=False, figsize=(25, 15))
# HGT_value = ko_B_percent_df.loc['HGT_pcofg'].tolist()
# print(HGT_value)
# print(ko_id_B_all)
# print(len(HGT_value))
# print(len(ko_id_B_all))
# marker = [u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o',u'o']
# color_list = ['r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r','r']
# color_list = 'rrrrrrrrrrrrrrrrbbbrrrrrbbbrrrrrrrrrrggggrrrrrrrrrrr'
# print(len(color_list))
# print(HGT_value)
# plt.plot(HGT_value, linewidth=0, marker='o', markersize=7, markeredgewidth=0, color=color_list)
# plt.plot(HGT_value, 'ro')
# plt.scatter(HGT_value, ko_id_B_all)
# plt.xticks(rotation=270)
# plt.savefig('/Users/songweizhi/Desktop/test_B.png', dpi=300)
# #print(len(ko_B_percent_df_NA.values))
