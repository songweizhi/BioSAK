import os
import glob
import argparse
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


annotation_file_dir = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber'
op_dir_filtered     = '/Users/songweizhi/PycharmProjects/BioSAK/demo_data/enrich/MAGs_dRep99_KEGG_D_GeneNumber_50'


annotation_file_re = '%s/*.txt' % annotation_file_dir
annotation_file_list = glob.glob(annotation_file_re)

interested_ko_set = set()
for each in open('/Users/songweizhi/Desktop/aaa.txt'):
    interested_ko_set.add(each.strip())
print(interested_ko_set)


for annotation_file in annotation_file_list:
    f_path, MAG_id, f_ext = sep_path_basename_ext(annotation_file)

    op_file = '%s/%s%s' % (op_dir_filtered, MAG_id.replace('_ko_stats_D_GeneNumber', ''), f_ext)

    op_file_handle = open(op_file, 'w')

    for each_line in open(annotation_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('KO	GeneNumber	Description'):
            op_file_handle.write(each_line)

        elif each_line_split[0] in interested_ko_set:
            op_file_handle.write(each_line)

    op_file_handle.close()





