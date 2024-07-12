import os
import argparse


srun_usage = '''
=============== srun example commands ===============

BioSAK srun -c "cd clock2_ns100000_run1; mcmctree"
BioSAK srun -c "cd clock2_ns100000_run2; mcmctree"

=====================================================
'''


def srun(args):

    cmd_to_run = args['c']

    cmd_to_run_split = cmd_to_run.strip().split(';')

    for each_cmd in cmd_to_run_split:
        each_cmd = each_cmd.strip()
        if each_cmd.startswith('cd '):
            os.chdir(each_cmd[3:])
        else:
            os.system(each_cmd)


if __name__ == '__main__':

    srun_parser = argparse.ArgumentParser()
    srun_parser.add_argument('-c', required=True,  help='command to run')
    args = vars(srun_parser.parse_args())
    srun(args)
