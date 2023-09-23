
#################################################### the simple one ####################################################

import argparse


parser = argparse.ArgumentParser()

parser.add_argument('-t', required=True, help='sequence file')
parser.add_argument('-i', required=False, type=int, help='mutation level')
parser.add_argument('-keep_cds', action="store_true", required=False, help='insert transfers only to non-coding regions')

args = vars(parser.parse_args())
input_seq_file_name = args['t']
mutation_level = args['i']
keep_cds = int(args['keep_cds'])


########################################### defined as required or optional ############################################

import argparse


parser = argparse.ArgumentParser(description='', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument('-h', action='help', help='Show this help message and exit')
required.add_argument('-snv_qc', dest='SNV_QC', nargs='?', required=True,  type=str, help='deepSNV QC file')
required.add_argument('-min_both', dest='MIN_BOTH', nargs='?', required=True, type=int, help='The minimum number of reads harboring SNV')

args = vars(parser.parse_args())
SNV_quality_file = args['SNV_QC']
min_var_reads_num = args['MIN_BOTH']


########################################### defined as required or optional ############################################

import argparse


parser = argparse.ArgumentParser(add_help=False)
subparsers = parser.add_subparsers(help="--", dest='subparser_name')


PrepIn_parser = subparsers.add_parser('PrepIn', formatter_class=CustomHelpFormatter, description='Runs tree, lineage_set, analyze, qa', epilog='Example: checkm lineage_wf ./bins ./output')
PrepIn_parser.add_argument('bin_folder', help="folder containing bins (fasta format)")
PrepIn_parser.add_argument('out_folder', help="folder to write output files")



