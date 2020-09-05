import os
import shutil
import argparse
import multiprocessing as mp


sra_reads_downloader_usage = '''
============================ sra_reads_downloader example commands =============================

module load sratoolkit/2.9.6-1
BioSAK sra_reads_downloader -in prokaryotes.csv

================================================================================================
'''


def dwnld_sra_reads_worker():

    pass


def dwnld_sra_reads(args):

    pass


args = {'i': 'SRR_Acc_List_1.txt',
        'o': 'test_out',
        't': 4}


# if __name__ == '__main__':
#
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('-in', required=True, help='input csv file')
#
#     args = vars(parser.parse_args())
#
#     sra_reads_downloader(args)
#
