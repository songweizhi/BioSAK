import os
import glob
import argparse
from Bio import SeqIO


cat_fa_usage = '''
================= cat_fa example commands =================

# example command 
BioSAK cat_fa -i seq_dir -x fa -o combined_seqs.fa

# combine sequence files in a folder after prefixing
# sequence ids with corresponding file name.

===========================================================
'''


def cat_fa(args):

    seq_folder    = args['i']
    file_ext      = args['x']
    sep_str       = args['s']
    combined_file = args['o']

    seq_file_re = '%s/*.%s' % (seq_folder, file_ext)
    seq_file_list = [os.path.basename(file_name) for file_name in glob.glob(seq_file_re)]

    combined_file_handle = open(combined_file, 'w')
    for seq_file in seq_file_list:
        pwd_seq_file = '%s/%s' % (seq_folder, seq_file)
        file_name_no_ext = seq_file[:-(len(file_ext)+1)]
        for each_seq in SeqIO.parse(pwd_seq_file, 'fasta'):
            seq_name_new = '%s%s%s' % (file_name_no_ext, sep_str, each_seq.id)
            each_seq.id = seq_name_new
            SeqIO.write(each_seq, combined_file_handle, 'fasta')
    combined_file_handle.close()


if __name__ == "__main__":

    cat_fa_parser = argparse.ArgumentParser()
    cat_fa_parser.add_argument('-i', required=True,                     help='sequence folder')
    cat_fa_parser.add_argument('-x', required=False, default='fasta',   help='sequence file extension, default: fasta')
    cat_fa_parser.add_argument('-s', required=False, default='__',      help='separator, default: __')
    cat_fa_parser.add_argument('-o', required=True,                     help='combined output file')
    args = vars(cat_fa_parser.parse_args())
    cat_fa(args)
