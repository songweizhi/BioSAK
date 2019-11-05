import os
import sys
from Bio import SeqIO

usage = """
FastaSplitler_by_size.py - splits a large FASTA file into smaller subsets by defined size(MB).

Usage:
FastaSplitler_by_size.py [Input_file] [Size_in_MB]
"""

if len(sys.argv) != 3:
    print(usage)
    exit(1)

to_split = sys.argv[1]
size_cutoff_MB_str = sys.argv[2]

# Make sure that the file provided is readable.
try:
    open(to_split, 'r')
except IOError:
    print("Error! " + to_split + " does not exist, or is not readable!")
    exit(1)

# make sure the second argument is a positive value
try:
    float(size_cutoff_MB_str)
except ValueError:
    print("Error! Please provide a positive value!")
    exit(1)
if float(size_cutoff_MB_str) <= 0:
    print("Error! Please provide a positive value!")
    exit(1)


def split_with_size(input_list, size_cutoff):
    new_list = []
    temp_list = []
    sum = 0
    for each in input_list:
        if len(each.seq) >= size_cutoff:
            if temp_list == []:
                # append current element
                temp_list.append(each)
                new_list.append(temp_list)
                # reset temp_list
                temp_list = []
            else:
                # if temp_list != [], append temp_list first
                new_list.append(temp_list)
                # reset temp_list
                temp_list = []
                # then, append new element
                temp_list.append(each)
                new_list.append(temp_list)
                # reset temp_list and sum
                temp_list = []
                sum = 0
        else:
            sum += len(each.seq)
            if sum <= size_cutoff:
                temp_list.append(each)
                pass
            else:
                new_list.append(temp_list)
                temp_list = []
                temp_list.append(each)
                sum = len(each.seq)
    new_list.append(temp_list)
    return new_list


size_cutoff_MB = float(size_cutoff_MB_str)
size_cutoff_bp = 1024 ** 2 * size_cutoff_MB

path, file_name = os.path.split(to_split)
file_name_no_extension = file_name.split('.')[0]
if path == '':
    path = '.'

all_fasta_list = list(SeqIO.parse(to_split, 'fasta'))
new_list = split_with_size(all_fasta_list, size_cutoff_bp)

n = 1
for each_subset in new_list:
    subset = open('%s/%s_%s.fasta' % (path, file_name_no_extension, n), 'w')
    for each_contig in each_subset:
        SeqIO.write(each_contig, subset, 'fasta')
    subset.close()
    n += 1

