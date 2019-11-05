#!/usr/bin/env python3

# Copyright (C) 2017, Weizhi Song, Kerrin Steensen and Torsten Thomas.
# songwz03@gmail.com and t.thomas@unsw.edu.au

# HgtSIM is free software: you can redistribute it and/or modify it under the terms of the GNU Affero
# General Public License as published by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.

# HgtSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General
# Public License for more details.

# You should have received a copy of the GNU Affero General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.


import os
import shutil
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def get_genome_size(fasta_file):
    genome = SeqIO.parse(fasta_file, 'fasta')
    total_length = 0
    for each_contig in genome:
        total_length += len(each_contig.seq)
    return total_length


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
    fasta_out.write_header()
    fasta_out.write_record(seq_record)
    fasta_out.write_footer()


def simulate_reads(pwd_genome_file, read_number, read_length, insert_size, split):
    path, file_name = os.path.split(pwd_genome_file)
    genome_name, ext = os.path.splitext(file_name)

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
        current_read_r1_id = 'r%s_from_%s_%sth_bp_#0/1' % (n, genome_name, rdm_num)
        current_read_r2_id = 'r%s_from_%s_%sth_bp_#0/2' % (n, genome_name, rdm_num)

        current_read_r1_id = '%s_%sth_r1' % (genome_name, n)
        current_read_r2_id = '%s_%sth_r2' % (genome_name, n)

        if split == 1:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_r1_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_r2_handle)
        else:
            export_dna_record(current_fragment_r1, current_read_r1_id, '', output_combined_handle)
            export_dna_record(current_fragment_r2_reverse_complement, current_read_r2_id, '', output_combined_handle)
        n += 1
    if split == 1:
        output_r1_handle.close()
        output_r2_handle.close()
    else:
        output_combined_handle.close()


def Reads_simulator(args):

    rederence_genome = args['R']
    reads_number = args['n']
    reads_length = args['l']
    insert_size = args['i']
    split = args['split']

    simulate_reads(rederence_genome, reads_number, reads_length, insert_size, split)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-R', required=True, help='reference genomes')
    parser.add_argument('-n', required=True, type=int, help='reads number')
    parser.add_argument('-l', required=True, type=int, help='reads length')
    parser.add_argument('-i', required=True, type=int, help='insert size')
    parser.add_argument('-split', action="store_true", help='Export forward and reverse reads to separate files')

    args = vars(parser.parse_args())

    Reads_simulator(args)
