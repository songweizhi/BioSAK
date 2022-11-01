import argparse
from Bio import SeqIO


def rm_short_reads(args):

    reads_file        = args['i']
    min_read_len      = args['l']
    output_reads_file = args['o']
    in_fastq          = args['fq']

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    output_reads_file_handle = open(output_reads_file, 'w')
    for each_read in SeqIO.parse(reads_file, seq_in_format):
        if len(each_read.seq) >= min_read_len:
            if seq_in_format == 'fasta':
                output_reads_file_handle.write('>%s\n' % each_read.id)
                output_reads_file_handle.write('%s\n' % str(each_read.seq))
            elif seq_in_format == 'fastq':
                SeqIO.write(each_read, output_reads_file_handle, 'fastq')
    output_reads_file_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=True,                        help='input file')
    parser.add_argument('-l',   required=True, type=int,              help='read length cutoff (bp)')
    parser.add_argument('-o',   required=True,                        help='output file')
    parser.add_argument('-fq',  required=False, action="store_true",  help='in fastq format, default: fa')
    args = vars(parser.parse_args())
    rm_short_reads(args)
