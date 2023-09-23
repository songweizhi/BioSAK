
########################################################################################################################

import os
import multiprocessing as mp


def exe_cmds(cmd_file, num_threads):

    cmd_set = set()
    for each_cmd in open(cmd_file):
        cmd_set.add(each_cmd.strip())
    
    print('Running %s commands with %s cores' % (len(cmd_set), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, cmd_set)
    pool.close()
    pool.join()
    print('Done!')


########################################################################################################################

import os
import glob
from Bio import SeqIO
import multiprocessing as mp


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def gbk2faa_worker(arg_list):

    gbk_in  = arg_list[0]
    faa_out = arg_list[1]

    faa_out_handle = open(faa_out, 'w')
    for seq_record in SeqIO.parse(gbk_in, 'genbank'):
        for feature in seq_record.features:
            if feature.type == 'CDS':
                feature_locus_tag = feature.qualifiers['locus_tag'][0]
                feature_translation = feature.qualifiers['translation'][0]
                faa_out_handle.write('>%s\n' % feature_locus_tag)
                faa_out_handle.write('%s\n' % feature_translation)
    faa_out_handle.close()


def gbk2faa(gbk_dir, gbk_ext, faa_dir, faa_ext, num_threads):

    gbk_file_re   = '%s/*.%s' % (gbk_dir, gbk_ext)
    gbk_file_list = glob.glob(gbk_file_re)

    list_of_arg_list = []
    for each_gbk in gbk_file_list:
        gbk_path, gbk_base, gbk_ext = sep_path_basename_ext(each_gbk)
        pwd_faa_out = '%s/%s.%s' % (faa_dir, gbk_base, faa_ext)
        list_of_arg_list.append([each_gbk, pwd_faa_out])

    print('Processing %s files with %s cores' % (len(list_of_arg_list), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(gbk2faa_worker, list_of_arg_list)
    pool.close()
    pool.join()
    print('Done')


gbk_folder          = '/Users/songweizhi/Desktop/gbk_files'
gbk_file_ext        = 'gbk'
output_faa_folder   = '/Users/songweizhi/Desktop/faa_files'
faa_file_ext        = 'faa'
num_threads         = 6
gbk2faa(gbk_folder, gbk_file_ext, output_faa_folder, faa_file_ext, num_threads)

