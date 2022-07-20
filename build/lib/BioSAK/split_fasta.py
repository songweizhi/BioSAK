import os
import argparse
from Bio import SeqIO
from global_functions import sep_path_basename_ext

split_fasta_usage = '''
============ split_fasta example command ============

BioSAK split_fasta -i input.fa -o output_dir
BioSAK split_fasta -i input.fa -o output_dir -n 6

=====================================================
'''


def split_fasta(args):

    fasta_in         = args['i']
    per_file_seq_num = args['n']
    output_dir       = args['o']

    fasta_in_path, fasta_in_basename, fasta_in_ext = sep_path_basename_ext(fasta_in)

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


if __name__ == '__main__':

    split_fasta_parser = argparse.ArgumentParser(usage=split_fasta_usage)
    split_fasta_parser.add_argument('-i', required=True,                       help='input fasta file')
    split_fasta_parser.add_argument('-n', required=False, type=int, default=1, help='sequence per file, default: 1')
    split_fasta_parser.add_argument('-o', required=True,                       help='output dir')
    args = vars(split_fasta_parser.parse_args())
    split_fasta(args)
