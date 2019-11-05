import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq

def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def sep_reads_by_barcode(fastq_in, sample_to_barcode_file, direction):

    # get barcode to sample dict
    barcode_to_sample_dict = {}
    for barcode in open(sample_to_barcode_file):

        barcode_split = barcode.strip().split(',')
        sample_id = barcode_split[0]
        barcode_seq = barcode_split[1]

        if direction == 'F':
            barcode_to_sample_dict[barcode_seq] = sample_id

        if direction == 'R':
            barcode_seq_reverse_complement = str(Seq(barcode_seq).reverse_complement())
            barcode_to_sample_dict[barcode_seq_reverse_complement] = sample_id


    for each_read in SeqIO.parse(fastq_in, 'fastq'):
        fastq_in_basename = '.'.join(fastq_in.split('.')[:-1])
        each_read_seq = str(each_read.seq)
        each_read_seq_first_10bp = each_read_seq[:15]

        current_barcode = ''
        current_read_sample = ''
        for barcode in barcode_to_sample_dict:
            if barcode in each_read_seq_first_10bp:
                current_barcode = barcode
                current_read_sample = barcode_to_sample_dict[barcode]

        if (current_barcode != '') and (current_read_sample != ''):
            barcode_pos_start = each_read_seq.index(current_barcode)
            barcode_pos_end = barcode_pos_start + len(current_barcode)
            read_without_barcode = each_read_seq[barcode_pos_end:]
            each_read_qual_str = each_read.format("fastq").split('\n')[3]
            read_qual_str_without_barcode = each_read_qual_str[barcode_pos_end:]

            # write out reads without barcode to corresponding files
            pwd_output_file = '%s/%s_%s.fastq' % (output_folder, fastq_in_basename, current_read_sample)
            output_file_handle = open(pwd_output_file, 'a')
            output_file_handle.write('@%s\n' % each_read.description)
            output_file_handle.write('%s\n' % read_without_barcode)
            output_file_handle.write('+\n')
            output_file_handle.write('%s\n' % read_qual_str_without_barcode)
            output_file_handle.close()


# define input and output files
wd = '/Users/songweizhi/Desktop/1-1_raw_data'
sample_to_barcode_file =    'sample_to_barcode.txt'
r1_fastq =                  'R1_total.fastq'
r2_fastq =                  'R2_total.fastq'
output_folder =             'separated_reads'


# go to working directory
os.chdir(wd)

# create output folder
force_create_folder(output_folder)

# separate reads by barcode
#sep_reads_by_barcode(r1_fastq, sample_to_barcode_file, 'F')
sep_reads_by_barcode(r2_fastq, sample_to_barcode_file, 'R')
