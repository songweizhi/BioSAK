import os
import glob
import shutil
import argparse
import multiprocessing as mp
from datetime import datetime


Prokka_parser_usage = '''
==================== Prokka example commands ====================

# for completed genome
BioSAK Prokka -i gnm_folder -x fa -t 6
BioSAK Prokka -i MAG_folder -x fa -t 6 -meta

==========================================================================
'''


def Prokka(argument_list):

    input_genome = argument_list[0]
    input_genome_folder = argument_list[1]
    meta_mode = argument_list[2]

    # # prepare command (according to Prokka)
    # input_genome_basename, input_genome_ext = os.path.splitext(input_genome)
    # pwd_input_genome = '%s/%s' % (input_genome_folder, input_genome)
    #
    # prodigal_cmd_meta = 'prodigal -f sco -q -c -m -g 11 -p meta -i %s -o %s' % (pwd_input_genome, pwd_output_sco)
    # prodigal_cmd_nonmeta = 'prodigal -f sco -q -c -m -g 11 -i %s -o %s' % (pwd_input_genome, pwd_output_sco)
    #
    # if meta_mode is True:
    #     prodigal_cmd = prodigal_cmd_meta
    # else:
    #     prodigal_cmd = prodigal_cmd_nonmeta
    #
    # os.system(prodigal_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',          required=True,                          help='input genome folder')
    parser.add_argument('-x',          required=True, default='fasta',         help='file extension')
    parser.add_argument('-d',          required=True,                          help='genome domain, Bacteria or Archaea')
    parser.add_argument('-meta',       required=False, action="store_true",    help='annotation mode for metagenome assembled genomes (MAGs)')
    parser.add_argument('-t',          required=False, type=int, default=1,    help='number of threads')
    args = vars(parser.parse_args())
    Prokka(args)


'''
prokka --compliant --cpus 12 --kingdom Bacteria --prefix contig_2560 --locustag contig_2560 --outdir contig_2560_prokka_wd contig_2560.fa
prokka --compliant --cpus 12 --kingdom Archaea --prefix NorthSea_bin022 --locustag NorthSea_bin022 --outdir NorthSea_bin022_prokka_wd NorthSea_bin022.fasta
--metagenome
'''