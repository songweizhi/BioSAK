import os
import glob
import argparse
from Bio import SeqIO


rename_seq_usage = '''
======================= rename_seqs.py example commands =======================

python3 rename_seqs.py -i gnm_dir -x fa -o op_dir

===============================================================================
'''

def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def rename_seqs(args):

    seq_file_in     = args['i']
    file_extension  = args['x']
    op_dir          = args['o']


    seq_in_re = '%s/*.%s' % (seq_file_in, file_extension)
    seq_in_list = glob.glob(seq_in_re)

    if len(seq_in_list) == 0:
        print('No sequence file detected, program exited!')
        exit()

    if os.path.isdir(op_dir) is True:
        print('Output folder detected, program exited: %s' % op_dir)
        exit()
    else:
        os.mkdir(op_dir)


    file_index = 1
    for each_file in seq_in_list:
        f_path, f_base, f_ext = sep_path_basename_ext(each_file)

        print('Processing %s/%s: %s' % (file_index, len(seq_in_list), f_base))

        file_out = '%s/%s%s' % (op_dir, f_base, f_ext)

        file_out_handle = open(file_out, 'w')
        seq_index = 1
        for each_seq in SeqIO.parse(each_file, 'fasta'):
            file_out_handle.write('>%s_%s\n' % (f_base, seq_index))
            file_out_handle.write('%s\n' % each_seq.seq)
            seq_index += 1
        file_out_handle.close()
        file_index += 1

    print('Done!')


if __name__ == '__main__':

    rename_seq_parser = argparse.ArgumentParser(usage=rename_seq_usage)
    rename_seq_parser.add_argument('-i',         required=True,                         help='input sequence folder')
    rename_seq_parser.add_argument('-x',         required=False, default='fasta',       help='file extension, default: fasta')
    rename_seq_parser.add_argument('-o',         required=True,                         help='output folder')
    args = vars(rename_seq_parser.parse_args())
    rename_seqs(args)
