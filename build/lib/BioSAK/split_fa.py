import os
import math
import argparse
from Bio import SeqIO

split_fa_usage = '''
============= split_fa example command =============

BioSAK split_fa -i input.fa -o output_dir
BioSAK split_fa -i input.fa -o output_dir -ns 200
BioSAK split_fa -i input.fa -o output_dir -nf 10

====================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def split_fa(args):

    fasta_in         = args['i']
    per_file_seq_num = args['ns']
    num_of_file      = args['nf']
    output_dir       = args['o']

    fasta_in_path, fasta_in_basename, fasta_in_ext = sep_path_basename_ext(fasta_in)

    if (per_file_seq_num is None) and (num_of_file is None):
        print('Please either use the "-ns" flag to define the number of sequence per file or the "-nf" flag to specify the number of files to be generated.')
        exit()
    elif (per_file_seq_num is not None) and (num_of_file is not None):
        print('Options "-ns" and "-nf" are incompatible, only specify one of them.')
        exit()
    elif (per_file_seq_num is not None) and (num_of_file is None):

        if os.path.isdir(output_dir) is True:
            print('Output folder detected, please provide a different name!')
            exit()
        else:
            os.mkdir(output_dir)

        n = 1
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            file_index = (n - 1) // per_file_seq_num + 1
            pwd_sub_file = '%s/%s_%s.fa' % (output_dir, fasta_in_basename, file_index)
            if per_file_seq_num == 1:
                pwd_sub_file = '%s/%s.fa' % (output_dir, seq_record.id)
            # write out sequence
            with open(pwd_sub_file, 'a') as pwd_sub_file_handle:
                pwd_sub_file_handle.write('>%s\n' % seq_record.id)
                pwd_sub_file_handle.write('%s\n'  % str(seq_record.seq))
            n += 1
    elif (per_file_seq_num is None) and (num_of_file is not None):

        if os.path.isdir(output_dir) is True:
            print('Output folder detected, please provide a different name!')
            exit()
        else:
            os.mkdir(output_dir)

        # get total number of sequences
        total_seq_num = 0
        for each_seq in SeqIO.parse(fasta_in, 'fasta'):
            total_seq_num += 1
        seq_num_per_file = math.ceil(total_seq_num/num_of_file)

        # write out
        seq_index = 1
        sub_file_seq_dict = dict()
        for seq_record in SeqIO.parse(fasta_in, 'fasta'):
            file_index = math.ceil(seq_index/seq_num_per_file)
            pwd_sub_file = '%s/%s_%s.fa' % (output_dir, fasta_in_basename, file_index)
            if pwd_sub_file not in sub_file_seq_dict:
                sub_file_seq_dict[pwd_sub_file] = dict()
            sub_file_seq_dict[pwd_sub_file][seq_record.id] = str(seq_record.seq)
            # write out sequence
            # with open(pwd_sub_file, 'a') as pwd_sub_file_handle:
            #     pwd_sub_file_handle.write('>%s\n' % seq_record.id)
            #     pwd_sub_file_handle.write('%s\n'  % str(seq_record.seq))
            seq_index += 1

        # write out sequence
        for each_file in sub_file_seq_dict:
            current_file_dict = sub_file_seq_dict[each_file]
            with open(each_file, 'w') as each_file_handle:
                for each_seq in sorted(list(current_file_dict.keys())):
                    each_file_handle.write('>%s\n' % each_seq)
                    each_file_handle.write('%s\n' % current_file_dict[each_seq])


if __name__ == '__main__':

    split_fa_parser = argparse.ArgumentParser(usage=split_fa_usage)
    split_fa_parser.add_argument('-i',  required=True,                          help='input fasta file')
    split_fa_parser.add_argument('-o',  required=True,                          help='output dir')
    split_fa_parser.add_argument('-ns', required=False, default=None, type=int, help='number of sequences per file')
    split_fa_parser.add_argument('-nf', required=False, default=None, type=int, help='number of files to be generated')
    args = vars(split_fa_parser.parse_args())
    split_fa(args)
