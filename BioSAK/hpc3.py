import os
import argparse


hpc3_usage = '''
========================= hpc3 example commands =========================

BioSAK hpc3 -q cpu-share -n iqtree_1 -c "iqtree -h"
BioSAK hpc3 -q cpu-share -t 1 -n iqtree_3 -c "iqtree -h"
BioSAK hpc3 -q cpu-share -t 12 -tpc 3 -n mcmctree -c mcmctree_cmds.txt
BioSAK hpc3 -q cpu -a oces -conda mybase -n iqtree_2 -c "iqtree -h"

# To use srun, you commands must NOT contain the double quote symbol (").

=========================================================================
'''


def hpc3(args):

    cmd                 = args['c']
    job_name            = args['n']
    email_address       = args['m']
    walltime            = args['wt']
    node_num            = args['node']
    core_num            = args['t']
    core_num_per_cmd    = args['tpc']
    queue_name          = args['q']
    conda_env           = args['conda']
    setting_a           = args['a']
    use_srun            = args['srun']

    js_file = '%s.sh' % job_name
    js_file_handle = open(js_file, 'w')
    js_file_handle.write('#!/bin/bash\n')

    if email_address is not None:
        js_file_handle.write('#SBATCH --mail-user=%s\n' % email_address)
        js_file_handle.write('#SBATCH --mail-type=begin\n')
        js_file_handle.write('#SBATCH --mail-type=end\n')

    js_file_handle.write('#SBATCH -t %s\n'          % walltime)
    js_file_handle.write('#SBATCH -N %s -n %s\n'    % (node_num, core_num))
    js_file_handle.write('#SBATCH -J %s\n'          % job_name)
    if queue_name in ['cpu', 'gpu', 'himem']:
        js_file_handle.write('#SBATCH -A %s\n'      % setting_a)
    js_file_handle.write('#SBATCH -p %s\n\n'        % queue_name)
    js_file_handle.write('cd $SLURM_SUBMIT_DIR\n\n')

    if conda_env is not None:
        js_file_handle.write('source ~/.bashrc\nconda activate %s\n'  % conda_env)

    if os.path.isfile(cmd) is False:
        js_file_handle.write('%s\n' % cmd)
    else:
        # write out to job script using srun
        cmd_list = []
        for each_cmd in open(cmd):
            if len(each_cmd.strip()) > 0:
                cmd_list.append(each_cmd.strip())

        for each_cmd in cmd_list[:-1]:
            if use_srun is False:
                js_file_handle.write('srun -n %s %s &\n' % (core_num_per_cmd, each_cmd))
            else:
                js_file_handle.write('srun -n %s BioSAK srun -c "%s" &\n' % (core_num_per_cmd, each_cmd))

        if use_srun is False:
            js_file_handle.write('srun -n %s %s\n' % (core_num_per_cmd, cmd_list[-1]))
        else:
            js_file_handle.write('srun -n %s BioSAK srun -c "%s"\n' % (core_num_per_cmd, cmd_list[-1]))

        js_file_handle.write('wait\n')

    js_file_handle.close()

    # submit jobscript
    os.system('sbatch %s' % js_file)


if __name__ == '__main__':

    hpc3_parser = argparse.ArgumentParser(usage=hpc3_usage)
    hpc3_parser.add_argument('-c',        required=True,                          help='command to submit')
    hpc3_parser.add_argument('-n',        required=True,                          help='job name')
    hpc3_parser.add_argument('-m',        required=False, default=None,           help='email address')
    hpc3_parser.add_argument('-wt',       required=False, default='23:59:59',     help='walltime, default: 23:59:59')
    hpc3_parser.add_argument('-node',     required=False, type=int, default=1,    help='number of node, default: 1')
    hpc3_parser.add_argument('-t',        required=False, type=int, default=12,   help='number of core, default: 12')
    hpc3_parser.add_argument('-tpc',      required=False, type=int, default=1,    help='number of core per command, default: 1')
    hpc3_parser.add_argument('-a',        required=False, default='boqianpy',     help='-A, boqianpy or oces, default: boqianpy')
    hpc3_parser.add_argument('-q',        required=True,                          help='queue, select from: cpu, cpu-share')
    hpc3_parser.add_argument('-conda',    required=False, default=None,           help='conda environment')
    hpc3_parser.add_argument('-srun',     required=False,  action='store_true',   help='use srun')
    args = vars(hpc3_parser.parse_args())
    hpc3(args)
