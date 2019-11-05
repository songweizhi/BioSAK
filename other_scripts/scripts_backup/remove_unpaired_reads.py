import os
import sys
from time import sleep
from Bio import SeqIO

usage = """
#####################################################################################

    remove_unpaired_reads.py - remove unpaired reads from 2 paired raw reads files.

    Usage:
    python3 remove_unpaired_reads.py [prefix_R1.fastq] [prefix_R2.fastq]

    Reminder: It needs about 3 hours to process 30GB input reads (two reads files
    in total, uncompressed).

#####################################################################################
"""

if len(sys.argv) != 3:
    print(usage)
    exit(1)

reads_1_file = sys.argv[1]
reads_2_file = sys.argv[2]


def get_fastq_reads_id_list(reads_file, file_extension):
    reads = SeqIO.parse(reads_file, file_extension)
    reads_id_list = []
    total_reads = 0
    for each_read in reads:
        reads_id_list.append(each_read.id)
        total_reads += 1
    return reads_id_list, total_reads


def get_new_reads_file(input_reads_file, common_reads_list, output_reads_file, file_extension):
    out = open(output_reads_file, 'w')
    input_reads = SeqIO.parse(input_reads_file, file_extension)
    for each_read in input_reads:
        if each_read.id in common_reads_list:
            SeqIO.write(each_read, out, file_extension)
    out.close()


# check file extension
path_file_1, file_1_name = os.path.split(reads_1_file)
path_file_2, file_2_name = os.path.split(reads_2_file)
file_1_name_split = file_1_name.split('.')
file_2_name_split = file_2_name.split('.')
file_1_only_name = file_1_name_split[0]
file_2_only_name = file_2_name_split[0]
file_1_extension = file_1_name_split[-1]
file_2_extension = file_2_name_split[-1]
if file_1_extension != file_2_extension:
    print('Please use same file extension for the 2 input reads files')
    exit()

# get output file name
if path_file_1 == '':
    path_file_1 = '.'
if path_file_2 == '':
    path_file_2 = '.'
file_extension = file_1_extension

file_1_out_name = '%s_only_paired.%s' % (file_1_only_name, file_extension)
file_2_out_name = '%s_only_paired.%s' % (file_2_only_name, file_extension)
pwd_file_1_out = '%s/%s' % (path_file_1, file_1_out_name)
pwd_file_2_out = '%s/%s' % (path_file_2, file_2_out_name)

# get first reads id list
print('Get reads id list from %s' % file_1_name)
sleep(0.5)
file1_reads_id_list, file1_total_reads = get_fastq_reads_id_list(reads_1_file, file_extension)
print('%s reads found in %s' % (file1_total_reads, file_1_name))
sleep(0.5)

# get second reads id list
print('Get reads id list from %s' % file_2_name)
sleep(0.5)
file2_reads_id_list, file2_total_reads = get_fastq_reads_id_list(reads_2_file, file_extension)
print('%s reads found in %s' % (file2_total_reads, file_2_name))
sleep(0.5)

# get common reads id list
print('Get paired reads list...')
sleep(0.5)
common_reads_list = set(file1_reads_id_list).intersection(file2_reads_id_list)
total_common_reads = len(common_reads_list)
print('%s paired reads found between %s and %s' % (total_common_reads, file_1_name, file_2_name))
sleep(0.5)

if (file1_total_reads == total_common_reads) and (file2_total_reads == total_common_reads):
    print('All reads in %s and %s are paired, program exited' % (file_1_name, file_2_name))
    exit()

# get output files
print('Remove unpaired reads from %s' % file_1_name)
sleep(0.5)
get_new_reads_file(reads_1_file, common_reads_list, pwd_file_1_out, file_extension)
print('Unpaired reads removed and saved paired reads to %s' % pwd_file_1_out)
sleep(0.5)
print('Remove unpaired reads from %s' % file_2_name)
sleep(0.5)
get_new_reads_file(reads_2_file, common_reads_list, pwd_file_2_out, file_extension)
print('Unpaired reads removed and saved paired reads to %s' % pwd_file_2_out)
sleep(0.5)
print('Done!')
