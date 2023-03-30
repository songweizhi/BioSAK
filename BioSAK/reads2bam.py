import os
import argparse
from BioSAK.global_functions import sep_path_basename_ext


reads2bam_usage = '''
=========================== reads2bam example commands ===========================

# Dependencies: bowtie2 and samtools

BioSAK reads2bam -ref ref.fa -r1 R1.fq -r2 R1.fq -fq -index -t 12
BioSAK reads2bam -ref ref.fa -r1 R1.fa -r2 R1.fa -up unpaired.fa -index -t 12
BioSAK reads2bam -ref ref.fa -up unpaired.fa -index -local -t 12 -tmp
BioSAK reads2bam -ref ref.fa -up unpaired_R1.fa,unpaired_R2.fa -index -t 12 -tmp

==================================================================================
'''


def reads2bam(args):

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

    ref_path, ref_basename, ref_ext = sep_path_basename_ext(ref_seq)
    cmd_bowtie2_build = 'bowtie2-build -f %s %s --threads %s' % (ref_seq, ref_basename, thread_num)

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
        cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -S %s.sam %s'          % (ref_basename, r1_seq, r2_seq, ref_basename, bowtie2_parameter)
    elif (r1_seq is not None) and (r2_seq is not None) and (unpaired_seq is not None):
        cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s.sam %s'    % (ref_basename, r1_seq, r2_seq, unpaired_seq, ref_basename, bowtie2_parameter)
    elif (r1_seq is None) and (r2_seq is None) and (unpaired_seq is not None):
        cmd_bowtie2 = 'bowtie2 -x %s -U %s -S %s.sam %s'                % (ref_basename, unpaired_seq, ref_basename, bowtie2_parameter)
    else:
        print('Please check your input reads files')
        exit()

    cmd_samtools_view  = 'samtools view -bS %s.sam -o %s.bam'    % (ref_basename, ref_basename)
    cmd_samtools_sort  = 'samtools sort %s.bam -o %s_sorted.bam' % (ref_basename, ref_basename)
    cmd_samtools_index = 'samtools index %s_sorted.bam'          % ref_basename

    if index_ref is True:
        os.system(cmd_bowtie2_build)
    os.system(cmd_bowtie2)
    os.system(cmd_samtools_view)
    os.system(cmd_samtools_sort)
    os.system(cmd_samtools_index)

    if keep_tmp is False:
        os.system('rm %s.sam' % ref_basename)
        os.system('rm %s.bam' % ref_basename)


if __name__ == '__main__':

    reads2bam_parser = argparse.ArgumentParser(usage=reads2bam_usage)
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

