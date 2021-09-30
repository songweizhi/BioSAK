#!/usr/bin/env python3

import os
import glob
import pandas
import argparse
import numpy as np

'''

cd /Users/songweizhi/Desktop
Python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_dbCAN.py -in MAGs_dRep99_dbCAN_stats_GeneNumber -in_percent -skip_1st_row
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i MAGs_dRep99_dbCAN_stats_GeneNumber_CAZy_cate.txt -o MAGs_dRep99_dbCAN_stats_GeneNumber_CAZy_cate.png

'''


def cazy_id_to_cate(cazy_id):

    cazy_cate = ''
    stop_now = False
    for each_element in cazy_id:
        if stop_now is False:
            if each_element.isalpha() is True:
                cazy_cate += each_element
            else:
                stop_now = True

    return cazy_cate


def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


def dict_to_file(genome_name_list, annotation_dict_cazy_cate, column_list_cazy_cate, in_percent, output_csv):
    # get category_num_lol
    category_num_lol = []
    for genome in genome_name_list:

        current_genome_annotation_results_dict = annotation_dict_cazy_cate[genome]

        current_genome_annotation_results_list = []
        for cog_cate in column_list_cazy_cate:
            cog_cate_num = 0
            if cog_cate in current_genome_annotation_results_dict:
                cog_cate_num = current_genome_annotation_results_dict[cog_cate]
            current_genome_annotation_results_list.append(cog_cate_num)

        # turn absolute number to percentage if specified
        if in_percent is True:
            current_genome_annotation_results_list_in_percent = turn_to_percentage(
                current_genome_annotation_results_list)
            category_num_lol.append(current_genome_annotation_results_list_in_percent)
        else:
            category_num_lol.append(current_genome_annotation_results_list)

    # turn list into arrary
    category_num_arrary = np.array(category_num_lol)

    # add row and column name to dataframe
    category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=column_list_cazy_cate)

    # write out
    category_num_df.to_csv(output_csv)


parser = argparse.ArgumentParser()
parser.add_argument('-in',            required=True,                          help='folder holds the func_stats.txt file for each genome')
parser.add_argument('-in_percent',    required=False, action="store_true",    help='in percent')
parser.add_argument('-skip_1st_row',  required=False, action="store_true",    help='skip 1st (header) row in input files')

args = vars(parser.parse_args())
input_folder =                  args['in']
in_percent =                    args['in_percent']
skip_1st_row =                  args['skip_1st_row']

if input_folder[-1] == '/':
    input_folder = input_folder[:-1]

output_csv_cazy_id   = '%s_CAZy_id.txt'   % input_folder
output_csv_cazy_cate = '%s_CAZy_cate.txt' % input_folder


input_files = '%s/*.txt' % input_folder
file_list = [os.path.basename(file_name) for file_name in glob.glob(input_files)]
file_list = sorted(file_list)


genome_name_list = []
detected_cazy_id_list = set()
detected_cazy_cate_list = set()
annotation_results_dict_cazy_id = {}
annotation_results_dict_cazy_cate = {}
cate_to_description_dict = {}
for each_file in file_list:
    genome_name = '_'.join(each_file.split('_')[:-2])
    genome_name_list.append(genome_name)
    current_annotation_results_dict = {}
    current_annotation_results_dict_cazy_cate = {}
    n = 0
    for each_category in open('%s/%s' % (input_folder, each_file)):
        each_category_split = each_category.strip().split('\t')
        cazy_id =           each_category_split[0]
        cazy_cate = cazy_id_to_cate(cazy_id)
        category_num_str =      each_category_split[1]
        category_description =  each_category_split[2]

        if skip_1st_row is True:
            if n > 0:
                current_annotation_results_dict[cazy_id] = int(category_num_str)
                detected_cazy_id_list.add(cazy_id)

                if cazy_cate not in current_annotation_results_dict_cazy_cate:
                    current_annotation_results_dict_cazy_cate[cazy_cate] = int(category_num_str)
                else:
                    current_annotation_results_dict_cazy_cate[cazy_cate] += int(category_num_str)
                detected_cazy_cate_list.add(cazy_cate)
        else:
            current_annotation_results_dict[cazy_id] = int(category_num_str)
            detected_cazy_id_list.add(cazy_id)

            if cazy_cate not in current_annotation_results_dict_cazy_cate:
                current_annotation_results_dict_cazy_cate[cazy_cate] = int(category_num_str)
            else:
                current_annotation_results_dict_cazy_cate[cazy_cate] += int(category_num_str)
            detected_cazy_cate_list.add(cazy_cate)

        if cazy_id not in cate_to_description_dict:
            cate_to_description_dict[cazy_id] = category_description
        n += 1

    annotation_results_dict_cazy_id[genome_name] = current_annotation_results_dict
    annotation_results_dict_cazy_cate[genome_name] = current_annotation_results_dict_cazy_cate

column_order_list_cazy_id = sorted([i for i in detected_cazy_id_list])
column_order_list_cazy_cate = sorted([i for i in detected_cazy_cate_list])

dict_to_file(genome_name_list, annotation_results_dict_cazy_id, column_order_list_cazy_id, in_percent, output_csv_cazy_id)
dict_to_file(genome_name_list, annotation_results_dict_cazy_cate, column_order_list_cazy_cate, in_percent, output_csv_cazy_cate)








