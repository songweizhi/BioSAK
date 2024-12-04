import os
import argparse
from distutils.spawn import find_executable


reads2bam_usage = '''
=========================== reads2bam example commands ===========================

# Dependencies: bowtie2 and samtools

BioSAK reads2bam -p sample1 -ref ref.fa -r1 R1.fq -r2 R1.fq -fq -index -t 12
BioSAK reads2bam -p sample2 -ref ref.fa -r1 R1.fa -r2 R1.fa -up unpaired.fa -index -t 12
BioSAK reads2bam -p sample3 -ref ref.fa -up unpaired.fa -index -local -t 12 -tmp
BioSAK reads2bam -p sample4 -ref ref.fa -up unpaired_R1.fa,unpaired_R2.fa -index -t 12 -tmp

==================================================================================
'''


def check_dependencies(program_list):
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)
    if not_detected_programs != []:
        print('%s not found, program exited!' % ','.join(not_detected_programs))
        exit()


def check_input_files (file_list):
    not_detected_files = []
    for each_input_file in file_list:
        if os.path.isfile(each_input_file) is False:
            not_detected_files.append(each_input_file)
    if not_detected_files != []:
        print('%s not found, program exited!' % ','.join(not_detected_files))
        exit()


def reads2bam(args):

    op_prefix       = args['p']
    ref_seq         = args['ref']
    index_ref       = args['index']
    r1_seq          = args['r1']
    r2_seq          = args['r2']
    unpaired_seq    = args['up']
    fq_format       = args['fq']
    local_aln       = args['local']
    no_unal         = args['no_unal']
    thread_num      = args['t']
    keep_tmp        = args['tmp']

    check_dependencies(['bowtie2-build', 'bowtie2', 'samtools'])
    check_input_files([ref_seq, r1_seq, r2_seq])

    cmd_bowtie2_build = 'bowtie2-build -f %s %s --threads %s --quiet' % (ref_seq, op_prefix, thread_num)
    bowtie2_parameter = ''
    if fq_format is True:
        bowtie2_parameter += ' -q'
    else:
        bowtie2_parameter += ' -f'

    if local_aln is True:
        bowtie2_parameter += ' --local'

    if no_unal is True:
        bowtie2_parameter += ' --no-unal'

    bowtie2_parameter += (' -p %s' % thread_num)

    cmd_bowtie2 = ''
    if (r1_seq is not None) and (r2_seq is not None) and (unpaired_seq is None):
        cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -S %s.sam %s'          % (op_prefix, r1_seq, r2_seq, op_prefix, bowtie2_parameter)
    elif (r1_seq is not None) and (r2_seq is not None) and (unpaired_seq is not None):
        cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s.sam %s'    % (op_prefix, r1_seq, r2_seq, unpaired_seq, op_prefix, bowtie2_parameter)
    elif (r1_seq is None) and (r2_seq is None) and (unpaired_seq is not None):
        cmd_bowtie2 = 'bowtie2 -x %s -U %s -S %s.sam %s'                % (op_prefix, unpaired_seq, op_prefix, bowtie2_parameter)
    else:
        print('Please check your input reads files')
        exit()

    cmd_samtools_view  = 'samtools view -bS %s.sam -o %s.bam'    % (op_prefix, op_prefix)
    cmd_samtools_sort  = 'samtools sort %s.bam -o %s_sorted.bam' % (op_prefix, op_prefix)
    cmd_samtools_index = 'samtools index %s_sorted.bam'          % op_prefix

    if index_ref is True:
        os.system(cmd_bowtie2_build)
    os.system(cmd_bowtie2)
    os.system(cmd_samtools_view)
    os.system(cmd_samtools_sort)
    os.system(cmd_samtools_index)

    if keep_tmp is False:
        os.system('rm %s.sam' % op_prefix)
        os.system('rm %s.bam' % op_prefix)


if __name__ == '__main__':

    reads2bam_parser = argparse.ArgumentParser(usage=reads2bam_usage)
    reads2bam_parser.add_argument('-p',         required=True,                          help='output prefix')
    reads2bam_parser.add_argument('-ref',       required=True,                          help='reference sequences')
    reads2bam_parser.add_argument('-index',     required=False, action="store_true",    help='index reference sequences')
    reads2bam_parser.add_argument('-r1',        required=False, default=None,           help='paired reads r1')
    reads2bam_parser.add_argument('-r2',        required=False, default=None,           help='paired reads r2')
    reads2bam_parser.add_argument('-up',        required=False, default=None,           help='unpaired reads')
    reads2bam_parser.add_argument('-fq',        required=False, action="store_true",    help='reads in fastq format')
    reads2bam_parser.add_argument('-local',     required=False, action="store_true",    help='perform local alignment')
    reads2bam_parser.add_argument('-no_unal',   required=False, action="store_true",    help='Suppress SAM records for reads that failed to align')
    reads2bam_parser.add_argument('-t',         required=False, type=int, default=1,    help='number of threads, default: 1')
    reads2bam_parser.add_argument('-tmp',       required=False, action="store_true",    help='keep temporary files')
    args = vars(reads2bam_parser.parse_args())
    reads2bam(args)

