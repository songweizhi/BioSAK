import os
import argparse


submitHPC_usage = '''
============= submitHPC example commands =============

BioSAK submitHPC -c cmds.txt -t 6 -p ale

# need to have submitHPC.sh in system path
submitHPC.sh --cmd "$cmd" -n 6 -c ale_1

======================================================
'''


def submitHPC(args):

    js_prefix   = args['p']
    cmds_txt    = args['c']
    thread_num  = args['t']

    js_index = 1
    for each_cmd in open(cmds_txt):
        cmd_str = each_cmd.strip()
        submit_cmd = 'submitHPC.sh --cmd "%s" -n %s -c %s_%s' % (cmd_str, thread_num, js_prefix, js_index)
        js_index += 1
        os.system(submit_cmd)


if __name__ == '__main__':

    # arguments for COG_parser
    submitHPC_parser = argparse.ArgumentParser(usage=submitHPC_usage)
    submitHPC_parser.add_argument('-p',               required=True,                              help='js prefix')
    submitHPC_parser.add_argument('-c',               required=True,                              help='commands txt')
    submitHPC_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads')
    args = vars(submitHPC_parser.parse_args())
    submitHPC(args)
