import os
import glob
import argparse
from Bio import SeqIO


gc_usage = '''
===== gc example commands ===== 

BioSAK gc -i mag_1.fa
BioSAK gc -i mag_dir -x fa 

===============================
'''


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def get_gc(sequence_file):

    total_gc_num = 0
    total_len = 0
    for ctg in SeqIO.parse(sequence_file, 'fasta'):
        seq_str = str(ctg.seq)
        gc_num = (seq_str.count('G') + seq_str.count('g') + seq_str.count('C') + seq_str.count('c'))
        total_gc_num += gc_num
        total_len += len(seq_str)

    gnm_gc = float("{0:.2f}".format(total_gc_num*100/total_len))

    return gnm_gc


def gc(argument_list):

    mag_folder =    argument_list['i']
    mag_extension = argument_list['x']

    if os.path.isfile(mag_folder) is True:
        _, f_base, _ = sep_path_basename_ext(mag_folder)
        gnm_gc = get_gc(mag_folder)
        print('%s\t%s' % (f_base, gnm_gc))

    elif os.path.isdir(mag_folder) is True:

        mag_file_re   = '%s/*.%s' % (mag_folder, mag_extension)
        mag_file_list = glob.glob(mag_file_re)

        if len(mag_file_list) == 0:
            print('No file detected, program exited')
            exit()

        for each_mag in sorted(mag_file_list):
            _, f_base, _ = sep_path_basename_ext(each_mag)
            gnm_gc = get_gc(each_mag)
            print('%s\t%s' % (f_base, gnm_gc))


if __name__ == '__main__':

    gc_parser = argparse.ArgumentParser()
    gc_parser.add_argument('-i',      required=True,                          help='sequence file/folder')
    gc_parser.add_argument('-x',      required=False,                         help='file extension')
    args = vars(gc_parser.parse_args())
    gc(args)
