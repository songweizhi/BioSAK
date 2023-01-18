import os
import glob
import argparse
from Bio import SeqIO
from itertools import groupby

from Bio import SeqIO


# total_len = 0
# for each_seq in SeqIO.parse('assembly_graph_NoEu.gfa_reads_min0bp.fastq', 'fastq'):
#     total_len += len(each_seq.seq)
#
# print(total_len)
# print(total_len/(1024*1024*1024))


def cigar_to_aln_len(cigar_string):
    # Given a CIGAR string, return the number of bases consumed from the query sequence
    read_consuming_ops = ("M", "I", "S", "=", "X")
    result = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in ("M", "=", "X"):
            print(op)
            result += length
    return result


cigar_to_aln_len('8H42M1I2M1D33M1D66M1D15M')