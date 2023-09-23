import os
import glob
import argparse


prefix_file_usage = '''
====================== prefix_file example commands ======================

BioSAK prefix_file -p NS -i refined_MAGs -x fa -o refined_MAGs_renamed

==========================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create) is True:
        os.system('rm -r %s' % folder_to_create)
    os.mkdir(folder_to_create)


def prefix_file(args):

    prefix          = args['p']
    dir_in          = args['i']
    file_ext        = args['x']
    dir_out         = args['o']
    force_overwrite = args['f']

    # create op dir
    if force_overwrite is True:
        if os.path.isdir(dir_out) is True:
            os.system('rm -r %s' % dir_out)
        os.system('mkdir %s' % dir_out)
    else:
        if os.path.isdir(dir_out) is False:
            os.system('mkdir %s' % dir_out)

    file_in_re   = '%s/*.%s' % (dir_in, file_ext)
    file_in_list = glob.glob(file_in_re)

    for each_file in file_in_list:
        f_path, f_base, f_ext   = sep_path_basename_ext(each_file)
        pwd_file_renamed        = '%s/%s_%s%s'  % (dir_out, prefix, f_base, f_ext)
        rename_cmd              = 'mv %s %s'    % (each_file, pwd_file_renamed)
        os.system(rename_cmd)


if __name__ == '__main__':

    prefix_file_parser = argparse.ArgumentParser()
    prefix_file_parser.add_argument('-p',   required=True,                          help='add prefix to sequence')
    prefix_file_parser.add_argument('-i',   required=True,                          help='input sequence file')
    prefix_file_parser.add_argument('-x',   required=True,                          help='file extension, default: fasta')
    prefix_file_parser.add_argument('-o',   required=True,                          help='the number of columns to keep')
    prefix_file_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(prefix_file_parser.parse_args())
    prefix_file(args)
