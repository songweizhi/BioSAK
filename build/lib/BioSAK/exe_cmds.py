import os
import argparse
import multiprocessing as mp

exe_cmds_usage = '''
====== exe_cmds example command ======

BioSAK exe_cmds -c cmds.txt -t 12

# one command per line
======================================
'''


def exe_cmds(args):

    cmd_file    = args['c']
    num_threads = args['t']

    cmd_set = set()
    for each_cmd in open(cmd_file):
        cmd_set.add(each_cmd.strip())

    print('Running %s commands with %s cores' % (len(cmd_set), num_threads))
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, sorted([i for i in cmd_set]))
    pool.close()
    pool.join()
    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c',     required=True,                       help='cmds file')
    parser.add_argument('-t',     required=False, type=int, default=1, help='number of threads')
    args = vars(parser.parse_args())
    exe_cmds(args)
