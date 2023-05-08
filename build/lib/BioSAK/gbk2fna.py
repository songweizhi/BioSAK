import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp


gbk2fna_usage = '''
========= gbk2fna example commands =========

BioSAK gbk2fna -i MAG.gbk -o MAG.fna
BioSAK gbk2fna -i gbk_dir -o fna_dir -t 6

============================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def gbk2fna_single(arg_list):
    gbk_in  = arg_list[0]
    fna_out = arg_list[1]
    SeqIO.convert(gbk_in, 'genbank', fna_out, 'fasta')


def gbk2fna(args):

    gbk_in      = args['i']
    in_ext      = args['ix']
    fna_out     = args['o']
    out_ext     = args['ox']
    num_threads = args['t']

    if os.path.isfile(gbk_in) is True:
        gbk2fna_single([gbk_in, fna_out])

    elif os.path.isdir(gbk_in) is True:

        gbk_file_re   = '%s/*.%s' % (gbk_in, in_ext)
        gbk_file_list = glob.glob(gbk_file_re)

        if len(gbk_file_list) == 0:
            print('No input gbk found, program exited!')
            exit()

        if os.path.isdir(fna_out) is True:
            print('output folder detected, program exited!')
            exit()
        os.system('mkdir %s' % fna_out)

        mp_arg_lol = []
        for each_gbk in gbk_file_list:
            gbk_path, gbk_base, gbk_ext = sep_path_basename_ext(each_gbk)
            pwd_fna_out = '%s/%s.%s' % (fna_out, gbk_base, out_ext)
            #gbk2fna_single([each_gbk, pwd_fna_out])
            mp_arg_lol.append([each_gbk, pwd_fna_out])

        pool = mp.Pool(processes=num_threads)
        pool.map(gbk2fna_single, mp_arg_lol)
        pool.close()
        pool.join()

    else:
        print('No input gbk found, program exited!')
        exit()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',  required=True,                   help='input gbk')
    parser.add_argument('-ix', required=False, default='gbk',   help='input file extension, default: gbk')
    parser.add_argument('-o',  required=True,                   help='output fna (contig sequences)')
    parser.add_argument('-ox', required=False, default='fna',   help='output file extension, default: fna')
    parser.add_argument('-t',  required=False, type=int, default=1,        help='number of threads')
    args = vars(parser.parse_args())
    gbk2fna(args)
