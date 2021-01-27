import os
import argparse


js_cmds_usage = '''
========================= js_cmds example commands =========================

BioSAK js_cmds -p Test -cmd cmds.txt -header js_header.txt -n 3 -force

============================================================================
'''


def js_cmds(args):

    js_prefix       = args['p']
    cmd_file        = args['cmds']
    auto_submit     = args['auto']
    js_header_file  = args['header']
    cmd_num_per_js  = args['n']
    js_folder_hpc   = args['js_hpc']
    force_overwrite = args['force']

    pwd_js_folder       = '%s/%s_js' % (os.getcwd(), js_prefix)

    if js_folder_hpc is None:
        js_folder_hpc = pwd_js_folder


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

    js_index = 1
    n = 1
    for each_cmd in open(cmd_file):

        pwd_js_previous = '%s/%s_%s.sh' % (pwd_js_folder, js_prefix, (js_index - 1))
        pwd_js          = '%s/%s_%s.sh' % (pwd_js_folder, js_prefix, js_index)
        pwd_js_hpc      = '%s/%s_%s.sh' % (js_folder_hpc, js_prefix, js_index)

        if os.path.isfile(pwd_js) is False:

            if js_index > 1:

                if auto_submit is True:
                    pwd_js_previous_handle = open(pwd_js_previous, 'a')
                    pwd_js_previous_handle.write('qsub %s\n' % pwd_js_hpc)
                    pwd_js_previous_handle.close()

            pwd_js_handle = open(pwd_js, 'w')
            pwd_js_handle.write(js_header_str)
            pwd_js_handle.write('\n')
            pwd_js_handle.write(each_cmd)
            pwd_js_handle.close()
        else:
            pwd_js_handle = open(pwd_js, 'a')
            pwd_js_handle.write(each_cmd)
            pwd_js_handle.close()

        if n % cmd_num_per_js == 0:
            js_index += 1
        n += 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-p',      required=True,                       help='js prefix')
    parser.add_argument('-cmds',   required=True,                       help='cmds file')
    parser.add_argument('-auto',   required=False, action="store_true", help='automatically submit next js')
    parser.add_argument('-header', required=True,                       help='js header')
    parser.add_argument('-n',      required=True, type=int,             help='number of cmds per js')
    parser.add_argument('-js_hpc', required=False, default=None,        help='Full path to js folder on HPC')
    parser.add_argument('-force',  required=False, action="store_true", help='force overwrite existing results')

    args = vars(parser.parse_args())
    js_cmds(args)
