#!/usr/bin/env python3

import sys
import warnings
import argparse


if __name__ == '__main__':

    ########################################################################################### initialize subparsers ############################################################################################

    # initialize the options parser
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')

    # Annotation modules
    COG2003_parser =                    subparsers.add_parser('COG2003',                    description='Wrapper for COG annotation (v2003)')
    COG2014_parser =                    subparsers.add_parser('COG2014',                    description='Wrapper for COG annotation (v2014)')

    ###################################################################################### define arguments for subparsers #######################################################################################

    COG2003_parser.add_argument('-i',                   required=True,                                  help='path to input sequences (in multi-fasta format)')
    COG2003_parser.add_argument('-x',                   required=False,                                 help='file extension')
    COG2003_parser.add_argument('-m',                   required=True,                                  help='The type of input sequences, "N" for "nucleotide", "P" for "protein"')
    COG2003_parser.add_argument('-db_dir',              required=True,                                  help='folder holds whog, fun.txt, cddid.tbl and Cog.* files')
    COG2003_parser.add_argument('-t',                   required=False, type=int, default=1,            help='number of threads')

    COG2014_parser.add_argument('-i',                   required=True,                                  help='path to input sequences (in multi-fasta format)')
    COG2014_parser.add_argument('-x',                   required=False,                                 help='file extension')
    COG2014_parser.add_argument('-m',                   required=True,                                  help='sequence type, "N/n" for "nucleotide", "P/p" for "protein"')
    COG2014_parser.add_argument('-depth',               required=False, default=None,                   help='gene depth file/folder')
    COG2014_parser.add_argument('-pct_by_all',          required=False, action='store_true',            help='normalize by all query genes, rather than those with COG assignment')
    COG2014_parser.add_argument('-db_dir',              required=True,                                  help='DB folder')
    COG2014_parser.add_argument('-diamond',             required=False, action='store_true',            help='run diamond (for big dataset), default is NCBI blastp')
    COG2014_parser.add_argument('-t',                   required=False, default=1,     type=int,        help='number of threads')
    COG2014_parser.add_argument('-evalue',              required=False, default=0.001, type=float,      help='evalue cutoff, default: 0.001')

    ###################################################### parse provided arguments and run corresponding function #####################################################

    # disable warning message
    warnings.filterwarnings('ignore')

    # parse options
    if (len(sys.argv) == 1) or (sys.argv[1] in ['-h', '-help', '--help']):
        print('haha')
        sys.exit(0)
    else:
        args = vars(parser.parse_args())

    if args['subparser_name'] == 'COG2003':
        print('COG2003')

    elif args['subparser_name'] == 'COG2014':
        print('COG2014')

    else:
        print('Unrecognized command: %s, program exited' % sys.argv[1])
        exit()


