import os
import argparse
from Bio import SeqIO

bam2reads_usage = '''
=========================== bam2reads example commands ===========================

module load samtools/1.15
BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg1 -r ctg_1
BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg2 -r ctg_2:200-5000
BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg3 -r ctg_3:200-
BioSAK bam2reads -b sorted.bam -o bam2reads_wd_ctg4 -r ctg_4:-500 -s reads.fastq

==================================================================================
'''


def select_seq(args):

    # read in argument
    seq_file      = args['seq']
    id_file       = args['id']
    select_option = args['option']
    output_file   = args['out']
    one_line      = args['oneline']
    in_fastq      = args['fq']

    # get provided id list
    seq_id_list = set()
    for seq_id in open(id_file):
        seq_id_list.add(seq_id.strip())

    seq_in_format = 'fasta'
    if in_fastq is True:
        seq_in_format = 'fastq'

    # extract sequences
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, seq_in_format):
        seq_id = seq_record.id
        if select_option == 1:
            if seq_id in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')
        if select_option == 0:
            if seq_id not in seq_id_list:

                if in_fastq is False:
                    if one_line is False:
                        SeqIO.write(seq_record, output_file_handle, 'fasta')
                    else:
                        SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
                else:
                    SeqIO.write(seq_record, output_file_handle, 'fastq')
    output_file_handle.close()


def bam2reads(args):

    bam_file        = args['b']
    region_str      = args['r']
    read_seq_file   = args['s']
    output_dir      = args['o']
    force_overwrite = args['force']

    # create output folder
    if os.path.isdir(output_dir) is True:
        if force_overwrite is False:
            print('Output folder detected, program exited!')
            exit()
        else:
            os.system('rm -r %s' % output_dir)
    os.system('mkdir %s' % output_dir)

    in_fastq = False
    if read_seq_file is not None:
        read_seq_file_ext = read_seq_file.split('.')[-1]
        if 'q' in read_seq_file_ext:
            in_fastq = True

    # file out
    mpileup_positions_file  = '%s/mpileup_positions.txt' % output_dir
    mpileup_op_file         = '%s/bam.mpileup'           % output_dir
    bam_header_txt          = '%s/bam_header.txt'        % output_dir
    read_id_txt             = '%s/reads_id.txt'          % output_dir
    extracted_read_file     = '%s/extracted_reads.fasta' % output_dir
    if in_fastq is True:
        extracted_read_file = '%s/extracted_reads.fastq' % output_dir

    # get bam header
    get_bam_header_cmd = 'samtools view -H %s > %s' % (bam_file, bam_header_txt)
    print('Running %s' % get_bam_header_cmd)
    os.system(get_bam_header_cmd)

    # get ref_seq_len_dict
    ref_seq_len_dict = {}
    for each_line in open(bam_header_txt):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            ref_seq_id = ''
            ref_seq_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    ref_seq_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_seq_len = int(each_element[3:])
            ref_seq_len_dict[ref_seq_id] = ref_seq_len

    # get region_list
    seq_id = region_str
    region_list = []
    if ':' in region_str:
        region_str_split = region_str.split(':')
        seq_id = region_str_split[0]
        region_list_original = region_str_split[1].split('-')
        ref_seq_len = ref_seq_len_dict[seq_id]

        # update left side pos
        if region_list_original[0] == '':
            region_list.append(1)
        else:
            region_list.append(int(region_list_original[0]))

        # update right side pos
        if region_list_original[1] == '':
            region_list.append(ref_seq_len)
        else:
            if int(region_list_original[1]) <= ref_seq_len:
                region_list.append(int(region_list_original[1]))
            else:
                region_list.append(ref_seq_len)
    else:
        region_list = [1, ref_seq_len_dict[seq_id]]

    # write out mpileup_positions_file
    mpileup_positions_file_handle = open(mpileup_positions_file, 'w')
    for each_pos in range(region_list[0], (region_list[1] + 1)):
        mpileup_positions_file_handle.write('%s\t%s\n' % (seq_id, each_pos))
    mpileup_positions_file_handle.close()

    # run samtools mpileup
    mpileup_cmd = 'samtools mpileup -Q 0 -l %s --output-QNAME -o %s %s' % (mpileup_positions_file, mpileup_op_file, bam_file)
    print('Running %s' % mpileup_cmd)
    os.system(mpileup_cmd)

    # parse samtools mpileup output
    print('Parsing mpileup output')
    aligned_read_set = set()
    for each_line in open(mpileup_op_file):
        each_line_split = each_line.strip().split('\t')
        current_aligned_read_list = each_line_split[6].split(',')
        for each_read in current_aligned_read_list:
            aligned_read_set.add(each_read)

    # write out read id
    read_id_txt_handle = open(read_id_txt, 'w')
    for each_read in aligned_read_set:
        read_id_txt_handle.write(each_read + '\n')
    read_id_txt_handle.close()
    print('ID of %s reads aligned to %s were exported to: %s' % (len(aligned_read_set), region_str, read_id_txt))

    # extract sequences
    if read_seq_file is not None:
        select_seq_arg_dict = {'seq': read_seq_file,  'id': read_id_txt,  'option': 1, 'out': extracted_read_file,  'oneline': True,  'fq': in_fastq}
        select_seq(select_seq_arg_dict)
    print('Sequence of %s reads aligned to %s were exported to: %s' % (len(aligned_read_set), region_str, extracted_read_file))
    print('Done!')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-b',     required=True,                       help='Sorted Bam file')
    parser.add_argument('-r',     required=True,                       help='Interested region')
    parser.add_argument('-s',     required=False, default=None,        help='Read sequence file')
    parser.add_argument('-o',     required=True,                       help='Output folder')
    parser.add_argument('-force', required=False, action="store_true", help='Force overwriting')
    args = vars(parser.parse_args())
    bam2reads(args)
