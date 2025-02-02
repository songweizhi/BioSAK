import os
import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


Reads_simulator_usage = '''
===================== Reads_simulator example commands =====================

# simulate 100000 pairs of reads from ref.fa
BioSAK Reads_simulator -p Test -r ref.fa -l 250 -i 300 -split -n 100000

# simulate reads to a depth of 50X
BioSAK Reads_simulator -p Test -r ref.fa -l 250 -i 300 -split -d 50

============================================================================
'''


def Reads_simulator(args):

    op_prefix       = args['p']
    pwd_genome_file = args['r']
    read_number     = args['n']
    read_depth      = args['d']
    read_length     = args['l']
    insert_size     = args['i']
    split           = args['split']

    if (read_number is None) and (read_depth is None):
        print('Please either provide the number of read pairs to simulate (-n) or the desired sequencing depth (-d).')
        print('Program exited!')
        exit()
    elif (read_number is not None) and (read_depth is not None):
        print('"-d" and "-n" can not be specified at the same time ')
        print('Program exited!')
        exit()

    # read in reference sequences
    ref_seq_total_len = 0
    ref_seq_dict = {}
    ref_seq_len_dict = {}
    for each_seq in SeqIO.parse(pwd_genome_file, 'fasta'):
        seq_id = each_seq.id
        seq_str = str(each_seq.seq)
        ref_seq_dict[seq_id] = seq_str
        ref_seq_len_dict[seq_id] = len(seq_str)
        ref_seq_total_len += len(seq_str)

    # calculate the number of reads to simulate
    if (read_number is None) and (read_depth is not None):
        read_number = round(ref_seq_total_len*read_depth/(read_length*2))
        print('The number of read pairs to simulate: %s' % read_number)

    # get the number of reads to simulate from each ref seq
    ref_to_read_num_dict = {}
    current_ref_index = 1
    total_assigned_read_num = 0
    for each_ref_seq in ref_seq_len_dict:
        if current_ref_index == len(ref_seq_len_dict):
            current_ref_seq_read_num = read_number - total_assigned_read_num
        else:
            current_ref_seq_len_pct = ref_seq_len_dict[each_ref_seq]/ref_seq_total_len
            current_ref_seq_read_num = round(current_ref_seq_len_pct*read_number)
            total_assigned_read_num += current_ref_seq_read_num
        ref_to_read_num_dict[each_ref_seq] = current_ref_seq_read_num
        current_ref_index += 1

    # simulate reads
    output_r1_handle = ''
    output_r2_handle = ''
    output_combined_handle = ''
    if split == 1:
        output_r1 = '%s_R1.fasta' % (op_prefix)
        output_r2 = '%s_R2.fasta' % (op_prefix)
        output_r1_handle = open(output_r1, 'w')
        output_r2_handle = open(output_r2, 'w')
    else:
        output_combined = '%s_R12.fasta' % (op_prefix)
        output_combined_handle = open(output_combined, 'w')

    overall_n = 1
    for each_ref_seq in ref_to_read_num_dict:

        ref_seq_str = ref_seq_dict[each_ref_seq]
        ref_seq_read_num = ref_to_read_num_dict[each_ref_seq]

        sequence_length = len(ref_seq_str)
        fragment_length = 2 * read_length + insert_size

        n = 1
        while n <= ref_seq_read_num:
            rdm_num = random.randint(1, sequence_length)
            current_fragment = ''
            if (rdm_num + fragment_length) <= sequence_length:
                current_fragment = ref_seq_str[rdm_num - 1: rdm_num + fragment_length - 1]
            elif (rdm_num + fragment_length) > sequence_length:
                seq_part_1_seq = ref_seq_str[rdm_num - 1:]
                seq_part_2_seq = ref_seq_str[:fragment_length - sequence_length + rdm_num - 1]
                current_fragment = seq_part_1_seq + seq_part_2_seq
            current_fragment_r1 = current_fragment[:read_length]
            current_fragment_r2 = current_fragment[-read_length:]
            current_fragment_r2_reverse_complement = str(Seq(current_fragment_r2).reverse_complement())
            current_read_r1_id = '%s_%s.1' % (op_prefix, overall_n)
            current_read_r2_id = '%s_%s.2' % (op_prefix, overall_n)

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
            overall_n += 1

    if split == 1:
        output_r1_handle.close()
        output_r2_handle.close()
        print('Simulated reads exported to %s and %s' % (output_r1, output_r2))
    else:
        output_combined_handle.close()
        print('Simulated reads exported to %s' % output_combined)

    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p',     required=True,                            help='output prefix')
    parser.add_argument('-r',     required=True,                            help='reference genome')
    parser.add_argument('-n',     required=False, type=int,   default=None, help='number of read pairs to simulate')
    parser.add_argument('-d',     required=False, type=float, default=None, help='desired read depth (X)')
    parser.add_argument('-l',     required=True,  type=int,                 help='read length (bp)')
    parser.add_argument('-i',     required=True,  type=int,                 help='insert size (bp)')
    parser.add_argument('-split', action="store_true",                      help='Export forward and reverse reads to separate files')
    args = vars(parser.parse_args())
    Reads_simulator(args)
