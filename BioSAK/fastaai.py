import os
import glob
import argparse
import subprocess


fastaai_usage = '''
============ fastaai example commands ============

BioSAK fastaai -i faa_files -x faa -o fastaai_op

==================================================
'''

def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


def fastaai(args):

    dir_in              = args['i']
    file_ext            = args['x']
    op_dir              = args['o']
    force_create_op_dir = args['f']
    num_threads         = args['t']

    check_executables(['fastaai'])

    faa_file_re   = '%s/*.%s' % (dir_in, file_ext)
    faa_file_list = glob.glob(faa_file_re)
    if len(faa_file_list) == 0:
        print('Sequence file not found, program exited!')
        exit()

    db_dir  = '%s/db'      % op_dir
    cmd_txt = '%s/cmd.txt' % op_dir

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output directory exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % db_dir)

    db_file      = '%s/database/user_genomes.db'                                                                             % db_dir
    cmd_build_db = 'fastaai build_db --proteins %s --verbose --output %s --database user_genomes.db --threads %s --compress' % (dir_in, db_dir, num_threads)
    cmd_fastaai  = 'fastaai db_query --query %s --target %s --verbose --output %s --threads %s --output_style matrix'        % (db_file, db_file, op_dir, num_threads)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'w')
    cmd_txt_handle.write(cmd_build_db + '\n')
    cmd_txt_handle.write(cmd_fastaai + '\n')
    cmd_txt_handle.close()
    print('Commands exported to %s' % cmd_txt)

    # run commands
    os.system(cmd_build_db)
    os.system(cmd_fastaai)

    print('Results exported to %s/results/FastAAI_matrix.txt' % op_dir)
    print('Done!')


if __name__ == '__main__':

    fastaai_parser = argparse.ArgumentParser(usage=fastaai_usage)
    fastaai_parser.add_argument('-i', required=True,                          help='faa files')
    fastaai_parser.add_argument('-x', required=False, default='faa',          help='file extension, default: faa')
    fastaai_parser.add_argument('-o', required=True,                          help='output directory')
    fastaai_parser.add_argument('-t', required=False, default=1, type=int,    help='number of threads, default: 1')
    fastaai_parser.add_argument('-f', required=False, action="store_true",    help='force overwrite')
    args = vars(fastaai_parser.parse_args())
    fastaai(args)
