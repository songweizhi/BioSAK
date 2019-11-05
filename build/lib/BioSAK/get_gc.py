import os
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_GC(gene_seq_file):

    seq_file_path, seq_file_basename, seq_file_extension = sep_path_basename_ext(gene_seq_file)
    output_file = '%s/%s_GC_content.txt' % (seq_file_path, seq_file_basename)

    output_file_handle = open(output_file, 'w')
    output_file_handle.write('ID\tGC(%)\n')
    for seq_record in SeqIO.parse(gene_seq_file, 'fasta'):

        seq_id = seq_record.id
        seq_sequence = str(seq_record.seq)
        seq_gc = (seq_sequence.count('G') + seq_sequence.count('g') + seq_sequence.count('C') + seq_sequence.count('c'))/len(seq_sequence)
        seq_gc = float("{0:.2f}".format(seq_gc*100))
        output_file_handle.write('%s\t%s\n' % (seq_id, seq_gc))

    output_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-in', required=True, help='input sequence file')

    args = vars(parser.parse_args())
    gene_seq_file = args['in']

    get_GC(gene_seq_file)
