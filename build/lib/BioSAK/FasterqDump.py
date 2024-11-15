import os
import glob
import argparse


FasterqDump_usage = '''
========== FasterqDump example commands ==========

BioSAK FasterqDump -i sra_dir -o op_dir -gzip -f

==================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def FasterqDump(args):

    sra_dir             = args['i']
    file_ext            = args['x']
    op_dir              = args['o']
    run_gzip            = args['gzip']
    force_create_op_dir = args['f']

    fasterq_dump_tmp    = '%s/tmp'                      % op_dir
    cmd_txt             = '%s/fasterq-dump_cmds.txt'    % op_dir

    sra_file_re = '%s/*.%s' % (sra_dir, file_ext)
    sra_file_list = glob.glob(sra_file_re)

    if len(sra_file_list) == 0:
        print('Input sra file not found, program exited!')
        exit()

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % fasterq_dump_tmp)

    fasterq_dump_cmd_list = []
    cmd_txt_handle = open(cmd_txt, 'w')
    for each_sra in sra_file_list:
        _, _, sra_base, _ = sep_path_basename_ext(each_sra)
        if run_gzip is False:
            fasterq_dump_cmd = 'fasterq-dump %s --split-3 -O %s -t %s' % (each_sra, sra_base, fasterq_dump_tmp)
        else:
            fasterq_dump_cmd = 'fasterq-dump %s --split-3 -O %s -t %s; gzip %s/%s_1.fastq; gzip %s/%s_2.fastq' % (each_sra, op_dir, fasterq_dump_tmp, op_dir, sra_base, op_dir, sra_base)
        cmd_txt_handle.write(fasterq_dump_cmd + '\n')
        fasterq_dump_cmd_list.append(fasterq_dump_cmd)
    cmd_txt_handle.close()

    # run fasterq-dump
    for each_cmd in fasterq_dump_cmd_list:
        print('Processing %s' % each_cmd.split(' ')[1])
        os.system(each_cmd)

    print('Done!')


if __name__ == '__main__':

    FasterqDump_parser = argparse.ArgumentParser()
    FasterqDump_parser.add_argument('-i',       required=True,                          help='sra directory')
    FasterqDump_parser.add_argument('-x',       required=False,default='sra',           help='sra directory, default is sra')
    FasterqDump_parser.add_argument('-o',       required=True,                          help='output directory')
    FasterqDump_parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    FasterqDump_parser.add_argument('-gzip',    required=False, action="store_true",    help='compress with gzip')
    args = vars(FasterqDump_parser.parse_args())
    FasterqDump(args)
