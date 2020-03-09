import os
import glob
import argparse
from BioSAK.global_functions import force_create_folder


split_folder_parser_usage = '''
================ split_folder example commands ================

BioSAK split_folder -in downloaded_genomes -x fasta -n 10

===============================================================
'''


def split_folder(args):

    file_folder = args['in']
    file_ext    = args['x']
    folder_num  = int(args['n'])

    # create folder
    file_folder_sep = '%s_sep' % file_folder
    force_create_folder(file_folder_sep)

    # get file list
    file_name_re = '%s/*.%s' % (file_folder, file_ext)
    file_name_list = [os.path.basename(file_name) for file_name in glob.glob(file_name_re)]
    file_name_list = sorted(file_name_list)

    # get the number of file per folder
    file_num_per_folder = round(len(file_name_list)/folder_num)

    n = 1
    while n <= folder_num:

        # define current folder name
        current_folder_name = '%s_%s' % (file_folder, n)

        # get file list in current folder
        if n < folder_num:
            files_in_current_folder = file_name_list[(file_num_per_folder * (n - 1)):(file_num_per_folder * n)]
        else:
            files_in_current_folder = file_name_list[(file_num_per_folder * (n - 1)):]

        # create folder
        pwd_current_folder = '%s/%s' % (file_folder_sep, current_folder_name)
        os.system('mkdir %s' % pwd_current_folder)

        # copy files to new folder
        for each_file in files_in_current_folder:
            os.system('cp %s/%s %s/' % (file_folder, each_file, pwd_current_folder))

        n += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('-in', required=True, help='file folder')
    parser.add_argument('-x', required=True, help='file extension')
    parser.add_argument('-n', required=False, type=int, help='number of subfolder')

    args = vars(parser.parse_args())

    split_folder(args)
