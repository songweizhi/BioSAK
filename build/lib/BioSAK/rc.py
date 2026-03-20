import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


rc_usage = '''
======== rc example commands ========

BioSAK rc -i AAAAATTTTTGGGGGCCCCC
BioSAK rc -i bin_1.ffn

=====================================
'''

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


def rc(args):

    seq_in = args['i']

    # check whether seq_in is a file
    if os.path.isfile(seq_in) is True:
        seq_in_path, seq_in_basename, seq_in_extension = sep_path_basename_ext(seq_in)
        seq_out = '%s/%s_rc%s' % (seq_in_path, seq_in_basename, seq_in_extension)
        seq_out_handle = open(seq_out, 'w')
        for seq_record in SeqIO.parse(seq_in, 'fasta'):
            seq_record_rc_id = '%s_rc' % seq_record.id
            seq_record_seq = seq_record.seq
            seq_record_seq_rc = seq_record_seq.reverse_complement()
            seq_out_handle.write('>%s\n' % seq_record_rc_id)
            seq_out_handle.write('%s\n' % str(seq_record_seq_rc))
        seq_out_handle.close()
    else:
        # check whether input is DNA sequence
        nc_bases = ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']
        nc_bases_count = 0
        for nc_base in nc_bases:
            nc_bases_count += seq_in.count(nc_base)

        if nc_bases_count == len(seq_in):
            seq_in_rc = str(Seq(seq_in).reverse_complement())
            print('>seq_in_rc\n%s\n' % seq_in_rc)
        else:
            print('Sequence file not exist or non DNA sequence detected, program exited!')
            exit()


if __name__ == '__main__':

    rc_parser = argparse.ArgumentParser(usage=rc_usage)
    rc_parser.add_argument('-i', required=True, help='input sequence(s)')
    args = vars(rc_parser.parse_args())
    rc(args)
