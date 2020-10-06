#!/usr/bin/env python3

import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


Reads_simulator_usage = '''
===================== Reads_simulator example commands =====================

To be added

============================================================================
'''


def get_genome_size(fasta_file):
    genome = SeqIO.parse(fasta_file, 'fasta')
    total_length = 0
    for each_contig in genome:
        total_length += len(each_contig.seq)
    return total_length


def Reads_simulator(args):

    pwd_genome_file = args['r']
    read_number     = args['n']
    read_length     = args['l']
    insert_size     = args['i']
    split           = args['split']

    path, file_name = os.path.split(pwd_genome_file)
    genome_name, ext = os.path.splitext(file_name)

    output_r1_handle = ''
    output_r2_handle = ''
    output_combined_handle = ''
    if split == 1:
        output_r1 = '%s_R1.fasta' % (genome_name)
        output_r2 = '%s_R2.fasta' % (genome_name)
        output_r1_handle = open(output_r1, 'w')
        output_r2_handle = open(output_r2, 'w')
    else:
        output_combined = '%s_R12.fasta' % (genome_name)
        output_combined_handle = open(output_combined, 'w')

    seq = str(SeqIO.read(pwd_genome_file, 'fasta').seq)
    sequence_length = len(seq)
    fragment_length = 2 * read_length + insert_size

    n = 1
    while n <= read_number:
        rdm_num = random.randint(1, sequence_length)
        current_fragment = ''
        if (rdm_num + fragment_length) <= sequence_length:
            current_fragment = seq[rdm_num - 1: rdm_num + fragment_length - 1]
        elif (rdm_num + fragment_length) > sequence_length:
            seq_part_1_seq = seq[rdm_num - 1:]
            seq_part_2_seq = seq[:fragment_length - sequence_length + rdm_num - 1]
            current_fragment = seq_part_1_seq + seq_part_2_seq
        current_fragment_r1 = current_fragment[:read_length]
        current_fragment_r2 = current_fragment[-read_length:]
        current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2, generic_dna).reverse_complement())
        current_read_r1_id = '%s_%s.1' % (genome_name, n)
        current_read_r2_id = '%s_%s.2' % (genome_name, n)

        if split == 1:
            output_r1_handle.write('>%s\n' % current_read_r1_id)
            output_r1_handle.write('%s\n'  % current_fragment_r1)
            output_r2_handle.write('>%s\n' % current_read_r2_id)
            output_r2_handle.write('%s\n'  % current_fragment_r2_reverse_complement)
        else:
            output_combined_handle.write('>%s\n' % current_read_r1_id)
            output_combined_handle.write('%s\n'  % current_fragment_r1)
            output_combined_handle.write('>%s\n' % current_read_r2_id)
            output_combined_handle.write('%s\n'  % current_fragment_r2_reverse_complement)

        n += 1

    if split == 1:
        output_r1_handle.close()
        output_r2_handle.close()
    else:
        output_combined_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-r',     required=True,           help='reference genomes')
    parser.add_argument('-n',     required=True, type=int, help='reads number')
    parser.add_argument('-l',     required=True, type=int, help='reads length')
    parser.add_argument('-i',     required=True, type=int, help='insert size')
    parser.add_argument('-split', action="store_true",     help='Export forward and reverse reads to separate files')

    args = vars(parser.parse_args())

    Reads_simulator(args)
