import os
import glob
import pandas
import argparse
import numpy as np


'''

Python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Kelp_dRep_fun_stats -out Kelp_dRep_fun_stats.txt -in_percent 
Python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in Tara_NM_fun_stats -out Tara_NM_fun_stats.txt -in_percent

Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i Kelp_dRep_fun_stats.txt -o Kelp_dRep_fun_stats.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i Tara_NM_fun_stats.txt -o Tara_NM_fun_stats.png

cd /Users/songweizhi/Desktop/validate_MetaCHIP
Python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in dRep95_fun_stat_files_sponge -out dRep95_fun_stat_sponge.txt -in_percent -skip_1st_row
Python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in dRep99_fun_stat_files_sponge -out dRep99_fun_stat_sponge.txt -in_percent -skip_1st_row

Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i dRep95_fun_stat_sponge.txt -o dRep95_fun_stat_sponge.png
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i dRep99_fun_stat_sponge.txt -o dRep99_fun_stat_sponge.png

'''


def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


parser = argparse.ArgumentParser()

parser.add_argument('-in',                          required=True,                          help='folder holds the func_stats.txt file for each genome')
parser.add_argument('-out',                         required=True,                          help='output csv file')
parser.add_argument('-columns',                     required=False,                         help='the order of columns, e.g. "J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S" for COG categories')
parser.add_argument('-in_percent',                  required=False, action="store_true",    help='in percent')
parser.add_argument('-skip_1st_row',                required=False, action="store_true",    help='skip 1st (header) row')
parser.add_argument('-with_functional_description', required=False, action="store_true",    help='add functional description to column header')

args = vars(parser.parse_args())
input_folder =                  args['in']
output_csv =                    args['out']
column_order =                  args['columns']
in_percent =                    args['in_percent']
skip_1st_row =                  args['skip_1st_row']
with_functional_description =   args['with_functional_description']

column_order_list = ['J', 'A', 'K', 'L', 'B', 'D', 'Y', 'V', 'T', 'M', 'N', 'Z', 'W', 'U', 'O', 'X', 'C', 'G', 'E', 'F', 'H', 'I', 'P', 'Q', 'R', 'S']

input_files = '%s/*.txt' % input_folder
file_list = [os.path.basename(file_name) for file_name in glob.glob(input_files)]
file_list = sorted(file_list)


############################################ get category_num_lol with dict ############################################

genome_name_list = []
detected_category_id_list = set()
annotation_results_dict = {}
cate_to_description_dict = {}
for each_file in file_list:
    genome_name = '_'.join(each_file.split('_')[:-2])
    genome_name_list.append(genome_name)

    current_annotation_results_dict = {}
    n = 0
    for each_category in open('%s/%s' % (input_folder, each_file)):
        each_category_split = each_category.strip().split('\t')
        category_id =           each_category_split[0]
        category_num_str =      each_category_split[1]
        category_description =  each_category_split[2]

        if skip_1st_row is True:
            if n > 0:
                current_annotation_results_dict[category_id] = int(category_num_str)
                detected_category_id_list.add(category_id)
        else:
            current_annotation_results_dict[category_id] = int(category_num_str)
            detected_category_id_list.add(category_id)

        if category_id not in cate_to_description_dict:
            cate_to_description_dict[category_id] = category_description

        n += 1

    annotation_results_dict[genome_name] = current_annotation_results_dict


# reorder columns

column_order_list_detected = []
for cate_id in column_order_list:
    if cate_id in detected_category_id_list:
        column_order_list_detected.append(cate_id)


column_order_list_detected_with_description = []
for detected_cate in column_order_list_detected:
    column_order_list_detected_with_description.append('%s_%s' % (detected_cate, cate_to_description_dict[detected_cate]))

# get category_num_lol
category_num_lol = []
for genome in genome_name_list:

    current_genome_annotation_results_dict = annotation_results_dict[genome]
    current_genome_annotation_results_list = []
    for cog_cate in column_order_list_detected:
        cog_cate_num = 0
        if cog_cate in current_genome_annotation_results_dict:
            cog_cate_num = current_genome_annotation_results_dict[cog_cate]
        current_genome_annotation_results_list.append(cog_cate_num)

    # turn absolute number to percentage if specified
    if in_percent is True:
        current_genome_annotation_results_list_in_percent = turn_to_percentage(current_genome_annotation_results_list)
        category_num_lol.append(current_genome_annotation_results_list_in_percent)
    else:
        category_num_lol.append(current_genome_annotation_results_list)


########################################################################################################################

# turn list into arrary
category_num_arrary = np.array(category_num_lol)

# add row and column name to dataframe
if with_functional_description is True:
    category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=column_order_list_detected_with_description)
else:
    category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=column_order_list_detected)

# write out
category_num_df.to_csv(output_csv)
