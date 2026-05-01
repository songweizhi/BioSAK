import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp


PebbleScout_usage = '''
================== PebbleScout example commands ==================

BioSAK PebbleScout -i gnm_dir -x fna -db meta,meta_vol2 -o op_dir

==================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def PebbleScout(args):

    file_dir            = args['i']
    file_ext            = args['x']
    db_str              = args['db']
    num_threads         = args['t']
    force_create_op_dir = args['f']
    op_dir              = args['o']

    # define file name
    gnm_dir             = '%s/genomes_concatenated' % op_dir
    pebblescout_op_dir  = '%s/pebblescout_outputs'  % op_dir
    cmd_txt             = '%s/cmds.txt'             % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gnm_dir)
    os.system('mkdir %s' % pebblescout_op_dir)

    db_list = db_str.split(',')

    # concatenate sequences in input genomes
    file_re = '%s/*.%s' % (file_dir, file_ext)
    file_list = glob.glob(file_re)
    for each_gnm in file_list:
        f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_gnm)
        current_gnm_after_cat = '%s/%s.%s' % (gnm_dir, f_base, f_ext)
        concatenates_seq_str = ''
        for each_seq in SeqIO.parse(each_gnm, 'fasta'):
            concatenates_seq_str = concatenates_seq_str + str(each_seq.seq).strip()
        current_gnm_after_cat_handle = open(current_gnm_after_cat, 'w')
        current_gnm_after_cat_handle.write('>%s\n' % f_base)
        current_gnm_after_cat_handle.write(concatenates_seq_str + '\n')
        current_gnm_after_cat_handle.close()

    # get pebblescout commands
    cmd_txt_handle = open(cmd_txt, 'w')
    cmd_set = set()
    for each_file in file_list:
        _, _, f_base, _ = sep_path_basename_ext(each_file)
        for each_db in db_list:
            pebblescout_cmd = 'curl -s -F "fasta=@%s" "https://pebblescout.ncbi.nlm.nih.gov/sra-cl-be/sra-cl-be.cgi?db=%s&m=2&rettype=pebblescout&download=yes" > %s/PebbleScout_%s_%s.txt' % (
                each_file, each_db, pebblescout_op_dir, f_base, each_db)
            cmd_set.add(pebblescout_cmd)
            cmd_txt_handle.write(pebblescout_cmd + '\n')
    cmd_txt_handle.close()

    # run pebblescout
    print('Running %s commands with %s cores' % (len(cmd_set), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, sorted([i for i in cmd_set]))
    pool.close()
    pool.join()
    print('Done!')


if __name__ == '__main__':

    PebbleScout_parser = argparse.ArgumentParser(usage=PebbleScout_usage)
    PebbleScout_parser.add_argument('-i',   required=True,                              help='input fasta file/dir')
    PebbleScout_parser.add_argument('-x',   required=False,default=None,                help='file extension')
    PebbleScout_parser.add_argument('-db',  required=False,default='meta,meta_vol2',    help='query database, default is meta and meta_vol2')
    PebbleScout_parser.add_argument('-t',   required=False, type=int, default=1,        help='number of core, default is 1')
    PebbleScout_parser.add_argument('-f',   required=False, action="store_true",        help='force overwrite')
    PebbleScout_parser.add_argument('-o',   required=True,                              help='output directory')
    args = vars(PebbleScout_parser.parse_args())
    PebbleScout(args)
