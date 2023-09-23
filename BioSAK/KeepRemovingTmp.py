import os
import time
import argparse
from datetime import datetime


KeepRemovingTmp_usage = '''
============== KeepRemovingTmp example command ==============

BioSAK KeepRemovingTmp -i gnm.fa -o ctg_id.txt
BioSAK KeepRemovingTmp -i reads.fastq -o reads_id.txt -fq

=============================================================
'''


def get_time_since_last_modification(target_folder_or_file):

    last_modified_time = os.path.getmtime(target_folder_or_file)
    current_time = time.time()
    tslm_sec  = current_time - last_modified_time
    tslm_min  = tslm_sec / 60
    tslm_hour = tslm_sec / (60 * 60)
    tslm_day  = tslm_sec / (60 * 60 * 24)
    tslm_sec  = float("{0:.2f}".format(tslm_sec))
    tslm_min  = float("{0:.2f}".format(tslm_min))
    tslm_hour = float("{0:.2f}".format(tslm_hour))
    tslm_day  = float("{0:.2f}".format(tslm_day))

    return tslm_sec, tslm_min, tslm_hour, tslm_day


def KeepRemovingTmp(args):

    tmp_dir  = args['d']
    time_str = args['l']

    time_format = '[%Y-%m-%d %H:%M:%S]'

    time_unit = time_str[-1]
    if time_unit not in ['s', 'm', 'h', 'd']:
        print('Invalid time format, program exited!')
        print('Here are a few examples: 30s, 30m, 6h, 5d')
        exit()

    time_value = float(time_str[:-1])

    time_len_sec = None
    if time_unit == 's':
        time_len_sec = time_value
    elif time_unit == 'm':
        time_len_sec = time_value*60
    elif time_unit == 'h':
        time_len_sec = time_value*60*60
    elif time_unit == 'd':
        time_len_sec = time_value*60*60*24

    sub_dir_list = next(os.walk(tmp_dir))[1]

    for each_sub_dir in sub_dir_list:
        pwd_sub_dir = '%s/%s' % (tmp_dir, each_sub_dir)
        subdir_tslm_sec, subdir_tslm_min, subdir_tslm_hour, subdir_tslm_day = get_time_since_last_modification(pwd_sub_dir)
        if subdir_tslm_sec >= time_len_sec:

            if time_unit == 's':
                print('%s Removing %s\t(last modified %ss ago)' % (datetime.now().strftime(time_format), each_sub_dir, subdir_tslm_sec))
            elif time_unit == 'm':
                print('%s Removing %s\t(last modified %sm ago)' % (datetime.now().strftime(time_format), each_sub_dir, subdir_tslm_min))
            elif time_unit == 'h':
                print('%s Removing %s\t(last modified %sh ago)' % (datetime.now().strftime(time_format), each_sub_dir, subdir_tslm_hour))
            elif time_unit == 'd':
                print('%s Removing %s\t(last modified %sd ago)' % (datetime.now().strftime(time_format), each_sub_dir, subdir_tslm_day))

            os.system('rm -r %s' % pwd_sub_dir)

    print('%s Done!' % datetime.now().strftime(time_format), )


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', required=True, help='tmp dir')
    parser.add_argument('-l', required=True, help='time length, e.g., 10s, 30m, 6h or 5d')
    args = vars(parser.parse_args())
    KeepRemovingTmp(args)


'''
python3 /home-user/wzsong/Scripts/rm_old_folder.py -d /home-user/wzsong/Japonicum/gapseq_metacyc/gapseq_tmp_dir -l 30m
'''
