import os
import math
import argparse


js_hpc3_usage = '''
==================== js_hpc3 example commands ====================

BioSAK js_hpc3 -p Demo -c cmds.txt -hd js_header.sh -n 3 -force

==================================================================
'''


def js_hpc3(args):

    js_prefix       = args['p']
    cmd_file        = args['c']
    js_header_file  = args['hd']
    cmd_num_per_js  = args['n']
    force_overwrite = args['f']
    cmd_t           = args['t']

    pwd_js_folder = '%s_job_scripts' % js_prefix

    if (os.path.isdir(pwd_js_folder) is True) and (force_overwrite is False):
        print('Output folder detected, program exited: %s' % pwd_js_folder)
        exit()
    else:
        if os.path.isdir(pwd_js_folder) is True:
            os.system('rm -r %s' % pwd_js_folder)
            os.mkdir(pwd_js_folder)
        else:
            os.mkdir(pwd_js_folder)

    # read in js header
    js_header_str = ''
    for each_line in open(js_header_file):
        js_header_str += each_line

    # get the total number of commands
    cmd_num = 0
    for each_cmd in open(cmd_file):
        if len(each_cmd.strip()) > 0:
            cmd_num += 1

    n = 1
    for each_cmd in open(cmd_file):
        if len(each_cmd.strip()) > 0:
            js_index = math.ceil((n-1)//cmd_num_per_js) + 1
            pwd_js = '%s/%s_%s.sh' % (pwd_js_folder, js_prefix, js_index)
            if os.path.isfile(pwd_js) is False:
                pwd_js_handle = open(pwd_js, 'w')
                pwd_js_handle.write(js_header_str)
                pwd_js_handle.write('\n')
                if n < cmd_num:
                    pwd_js_handle.write('srun -n %s %s &\n' % (cmd_t, each_cmd.strip()))
                else:
                    pwd_js_handle.write('%s\n' % each_cmd.strip())
                pwd_js_handle.close()
            else:
                pwd_js_handle = open(pwd_js, 'a')
                if (math.ceil((n-1)//cmd_num_per_js) + 1) == (math.ceil((n)//cmd_num_per_js) + 1):
                    if n < cmd_num:
                        pwd_js_handle.write('srun -n %s %s &\n' % (cmd_t, each_cmd.strip()))
                    else:
                        pwd_js_handle.write('srun -n %s %s\n' % (cmd_t, each_cmd.strip()))
                        pwd_js_handle.write('wait\n')
                else:
                    pwd_js_handle.write('srun -n %s %s\n' % (cmd_t, each_cmd.strip()))
                    pwd_js_handle.write('wait\n')
                pwd_js_handle.close()
            n += 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p',       required=True,                          help='js prefix')
    parser.add_argument('-c',       required=True,                          help='cmds file')
    parser.add_argument('-hd',      required=True,                          help='js header')
    parser.add_argument('-n',       required=True,  type=int,               help='number of cmds per job script')
    parser.add_argument('-t',       required=False, type=int, default=1,    help='number of threads specified in the commands')
    parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite existing results')
    args = vars(parser.parse_args())
    js_hpc3(args)
