import os
import glob
import pandas
import argparse
import numpy as np


'''

Python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_COG.py -in Kelp_dRep_fun_stats -out Kelp_dRep_fun_stats.txt -in_percent -columns J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S
Python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_COG.py -in Tara_NM_fun_stats -out Tara_NM_fun_stats.txt -in_percent -columns J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S

Rscript ~/PycharmProjects/MyBioTools/MyBioTools/COG_boxplot_last1row.R -i Kelp_dRep_fun_stats.txt -o Kelp_dRep_fun_stats.png
Rscript ~/PycharmProjects/MyBioTools/MyBioTools/COG_boxplot_last1row.R -i Tara_NM_fun_stats.txt -o Tara_NM_fun_stats.png

cd /Users/songweizhi/Desktop/555
python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_COG.py -in Kelp_and_HGT_COG_func_stats -out Kelp_and_HGT_COG_func_stats_num.txt -skip_1st_row -columns J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S
python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_COG.py -in Tara_NM_and_HGT_COG_func_stats -out Tara_NM_and_HGT_COG_func_stats_num.txt -skip_1st_row -columns J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S

python3 ~/PycharmProjects/MyBioTools/MyBioTools/boxplot_matrix_COG.py -in test -out test.txt -skip_1st_row -columns J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S

'''


def turn_to_percentage(number_list):
    number_list_percent = []
    for each_element in number_list:
        each_element_percent = float("{0:.2f}".format(each_element / sum(number_list)))
        number_list_percent.append(each_element_percent)
    return number_list_percent


parser = argparse.ArgumentParser()

parser.add_argument('-in',          required=True,                          help='folder holds the func_stats.txt file for each genome')
parser.add_argument('-out',         required=True,                          help='output csv file')
parser.add_argument('-columns',     required=False,                         help='the order of columns, e.g. "J,A,K,L,B,D,Y,V,T,M,N,Z,W,U,O,X,C,G,E,F,H,I,P,Q,R,S" for COG categories')
parser.add_argument('-in_percent',  required=False, action="store_true",    help='in percent')
parser.add_argument('-skip_1st_row',required=False, action="store_true",    help='skip 1st (header) row')

args = vars(parser.parse_args())
input_folder = args['in']
output_csv = args['out']
column_order = args['columns']
in_percent = args['in_percent']
skip_1st_row = args['skip_1st_row']


column_order_list = []
if column_order is not None:
    column_order_list = column_order.split(',')


input_files = '%s/*.txt' % input_folder
file_list = [os.path.basename(file_name) for file_name in glob.glob(input_files)]
file_list = sorted(file_list)


############################################ get category_num_lol with dict ############################################

genome_name_list = []
detected_category_id_list = set()
annotation_results_dict = {}
for each_file in file_list:
    genome_name = '_'.join(each_file.split('_')[:-2])
    genome_name_list.append(genome_name)

    current_annotation_results_dict = {}
    n = 0
    for each_category in open('%s/%s' % (input_folder, each_file)):
        each_category_split = each_category.strip().split('\t')
        category_id = each_category_split[0]
        category_num_str = each_category_split[1]

        if skip_1st_row is True:
            if n > 0:
                current_annotation_results_dict[category_id] = int(category_num_str)
                detected_category_id_list.add(category_id)
        else:
            current_annotation_results_dict[category_id] = int(category_num_str)
            detected_category_id_list.add(category_id)

        n += 1

    annotation_results_dict[genome_name] = current_annotation_results_dict


# reorder columns
if column_order is None:
    detected_category_id_list_reordered = sorted([i for i in detected_category_id_list])
else:
    detected_category_id_list_reordered = []
    for cate_id in column_order_list:
        if cate_id in detected_category_id_list:
            detected_category_id_list_reordered.append(cate_id)


# get category_num_lol
category_num_lol = []
for genome in genome_name_list:

    current_genome_annotation_results_dict = annotation_results_dict[genome]
    current_genome_annotation_results_list = []
    for cog_cate in detected_category_id_list_reordered:
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
category_num_df = pandas.DataFrame(category_num_arrary, index=genome_name_list, columns=detected_category_id_list_reordered)

# write out
category_num_df.to_csv(output_csv)
