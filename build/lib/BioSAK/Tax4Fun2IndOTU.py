import os
import glob
import argparse


Tax4Fun2IndOTU_usage = '''
============================== Tax4Fun2IndOTU example commands ==============================

=============================================================================================
'''


def Tax4Fun2IndOTU(args):

    print(args)


if __name__ == '__main__':

    COG_parser = argparse.ArgumentParser(usage=Tax4Fun2IndOTU_usage)
    COG_parser.add_argument('-i',               required=True,                              help='path to input sequences (in multi-fasta format)')
    COG_parser.add_argument('-x',               required=False,                             help='file extension')
    COG_parser.add_argument('-m',               required=True,                              help='sequence type, "N/n" for "nucleotide", "P/p" for "protein"')
    COG_parser.add_argument('-depth',           required=False, default=None,               help='gene depth file/folder')
    COG_parser.add_argument('-pct_by_all',      required=False, action='store_true',        help='normalize by all query genes, including those without COG assignment')
    COG_parser.add_argument('-db_dir',          required=True,                              help='DB folder')
    COG_parser.add_argument('-diamond',         required=False, action='store_true',        help='run diamond (for big dataset), default is NCBI blastp')
    COG_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads')
    COG_parser.add_argument('-evalue',          required=False, default=0.001, type=float,  help='evalue cutoff, default: 0.001')
    args = vars(COG_parser.parse_args())
    Tax4Fun2IndOTU(args)
