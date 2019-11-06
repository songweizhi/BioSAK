import os
import argparse
from time import sleep
from datetime import datetime


get_Pfam_hmms_usage = '''
===================================== get_Pfam_hmms example commands =====================================

# module needed
module load hmmer/3.2.1

# for completed genome
BioSAK get_Pfam_hmms -pfam Pfam-A.hmm -id needed_ids.txt

==========================================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_Pfam_hmms(args):

    downloaded_pfam_db = args['pfam']
    hmm_id_file = args['id']

    # define file name
    time_format = '[%Y-%m-%d %H:%M:%S]'
    pfam_db_path, pfam_db_basename, pfam_db_extension =             sep_path_basename_ext(downloaded_pfam_db)
    hmm_id_file_path, hmm_id_file_basename, hmm_id_file_extension = sep_path_basename_ext(hmm_id_file)
    downloaded_pfam_db_stats =  '%s/%s.stat'                        % (pfam_db_path, pfam_db_basename)
    Pfam_41_ids_updated =       '%s/%s_latest_version.txt'          % (hmm_id_file_path, hmm_id_file_basename)
    Pfam_profiles_updated =     '%s/%s.hmm'                         % (hmm_id_file_path, hmm_id_file_basename)

    # read in 41 Pfam ids
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Read in provided profile ids'))
    Pfam_41_id_list_no_version = []
    for each_hmm in open(hmm_id_file):
        Pfam_41_id_list_no_version.append(each_hmm.strip())

    # get all hmm profile id from downloaded db
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Get summary statistics for %s' % downloaded_pfam_db))
    hmmstat_cmd_db = 'hmmstat %s > %s' % (downloaded_pfam_db, downloaded_pfam_db_stats)
    os.system(hmmstat_cmd_db)

    # parse Pfam-A.hmm.stat and get updated id file
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Get the latest version of provided profiles'))
    id_to_version_dict_new = {}
    Pfam_41_ids_updated_handle = open(Pfam_41_ids_updated, 'w')
    for each_hmm in open(downloaded_pfam_db_stats):

        if not each_hmm.startswith('#') and (each_hmm.strip() != ''):

            each_hmm_split = each_hmm.strip().split(' ')

            # remove '' from each_hmm_split
            each_hmm_split_no_space = []
            for i in each_hmm_split:
                if i != '':
                    each_hmm_split_no_space.append(i)

            hmm_id_version = each_hmm_split_no_space[2]
            hmm_id_no_version = hmm_id_version.split('.')[0]
            id_to_version_dict_new[hmm_id_no_version] = hmm_id_version

            if hmm_id_no_version in Pfam_41_id_list_no_version:
                Pfam_41_ids_updated_handle.write('%s\n' % hmm_id_version)

    Pfam_41_ids_updated_handle.close()

    # extract updated hmm profile from downloaded db
    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Extract profiles'))
    hmmfetch_cmd = 'hmmfetch -f %s %s > %s' % (downloaded_pfam_db, Pfam_41_ids_updated, Pfam_profiles_updated)
    os.system(hmmfetch_cmd)

    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Delete temporary files'))
    os.system('rm %s' % downloaded_pfam_db_stats)

    sleep(0.5)
    print('%s %s' % ((datetime.now().strftime(time_format)), 'Done!'))


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for PI
    parser.add_argument('-pfam', required=True, help='Pfam db file, normally with name Pfam-A.hmm')
    parser.add_argument('-id',   required=True, help='ids of profiles need to be extracted, one id per line')

    args = vars(parser.parse_args())

    get_Pfam_hmms(args)
