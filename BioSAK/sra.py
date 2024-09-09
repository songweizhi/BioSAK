import os
import argparse
import multiprocessing as mp


sra_usage = '''
============= sra example commands =============

Requirement: sratoolkit

BioSAK sra -i sra_id.txt -o op_dir -t 12

================================================
'''

def sra(args):

    id_txt              = args['i']
    op_dir              = args['o']
    num_threads         = args['t']
    force_create_op_dir = args['f']

    fasterq_dump_tmp    = '%s/fasterq_dump_tmp' % op_dir
    cmd_txt             = '%s/cmds.txt'         % op_dir

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % fasterq_dump_tmp)

    id_set = set()
    for each_id in open(id_txt):
        id_set.add(each_id.strip().split()[0])

    prefetch_cmd_list = []
    fasterq_dump_cmd_list = []
    for each_id in id_set:
        prefetch_cmd     = 'prefetch %s -O %s'                            % (each_id, op_dir)
        fasterq_dump_cmd = 'fasterq-dump %s/%s --split-files -O %s -t %s' % (op_dir, each_id, op_dir, fasterq_dump_tmp)
        prefetch_cmd_list.append(prefetch_cmd)
        fasterq_dump_cmd_list.append(fasterq_dump_cmd)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'w')
    for each in sorted(prefetch_cmd_list):
        cmd_txt_handle.write(each + '\n')
    for each in sorted(fasterq_dump_cmd_list):
        cmd_txt_handle.write(each + '\n')
    cmd_txt_handle.close()

    # download files with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, prefetch_cmd_list)
    pool.close()
    pool.join()

    # extract fastq file from SRA with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, fasterq_dump_cmd_list)
    pool.close()
    pool.join()

    print('Done!')


if __name__ == '__main__':

    sra_parser = argparse.ArgumentParser()
    sra_parser.add_argument('-i',   required=True,                          help='id')
    sra_parser.add_argument('-o',   required=True,                          help='output directory')
    sra_parser.add_argument('-t',   required=False, type=int, default=1,    help='number of threads, default is 1')
    sra_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(sra_parser.parse_args())
    sra(args)
