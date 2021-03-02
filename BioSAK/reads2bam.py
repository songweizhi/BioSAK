import os
import argparse
from BioSAK.global_functions import sep_path_basename_ext


reads2bam_usage = '''
=================================== reads2bam example commands ===================================

module load bowtie/2.3.5.1
module load samtools/1.10
BioSAK reads2bam -p Demo -ref ref.fa -r1 R1.fa -r2 R1.fa -u unpaired.fa -index_ref -t 12
BioSAK reads2bam -p Demo -ref ref.fa -r1 R1.fq -r2 R1.fq -fastq -index_ref -t 12
BioSAK reads2bam -p Demo -ref ref.fa -u unpaired.fa -index_ref -t 12 -tmp

==================================================================================================
'''


def reads2bam(args):

    output_prefix   = args['p']
    ref_seq         = args['ref']
    index_ref       = args['index_ref']
    r1_seq          = args['r1']
    r2_seq          = args['r2']
    unpaired_seq    = args['u']
    fq_format       = args['fastq']
    thread_num      = args['t']
    keep_tmp        = args['tmp']


    ref_path, ref_basename, ref_ext = sep_path_basename_ext(ref_seq)

    cmd_bowtie2_build   = 'bowtie2-build -f %s %s --threads %s' % (ref_seq, ref_basename, thread_num)

    cmd_bowtie2 = ''
    if (r1_seq is not None) and (r2_seq is not None) and (unpaired_seq is None):
        cmd_bowtie2     = 'bowtie2 -x %s -1 %s -2 %s -S %s.sam -p %s -f' % (ref_basename, r1_seq, r2_seq, output_prefix, thread_num)
        if fq_format is True:
            cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -S %s.sam -p %s -q' % (ref_basename, r1_seq, r2_seq, output_prefix, thread_num)

    elif (r1_seq is not None) and (r2_seq is not None) and (unpaired_seq is not None):
        cmd_bowtie2     = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s.sam -p %s -f' % (ref_basename, r1_seq, r2_seq, unpaired_seq, output_prefix, thread_num)
        if fq_format is True:
            cmd_bowtie2 = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s.sam -p %s -q' % (ref_basename, r1_seq, r2_seq, unpaired_seq, output_prefix, thread_num)

    elif (r1_seq is None) and (r2_seq is None) and (unpaired_seq is not None):
        cmd_bowtie2     = 'bowtie2 -x %s -U %s -S %s.sam -p %s -f' % (ref_basename, unpaired_seq, output_prefix, thread_num)
        if fq_format is True:
            cmd_bowtie2 = 'bowtie2 -x %s -U %s -S %s.sam -p %s -q' % (ref_basename, unpaired_seq, output_prefix, thread_num)
    else:
        print('Please check your input reads files')
        exit()

    cmd_samtools_view   = 'samtools view -bS %s.sam -o %s.bam' % (output_prefix, output_prefix)
    cmd_samtools_sort   = 'samtools sort %s.bam -o %s_sorted.bam' % (output_prefix, output_prefix)
    cmd_samtools_index  = 'samtools index %s_sorted.bam' % output_prefix

    if index_ref is True:
        os.system(cmd_bowtie2_build)
    os.system(cmd_bowtie2)
    os.system(cmd_samtools_view)
    os.system(cmd_samtools_sort)
    os.system(cmd_samtools_index)

    if keep_tmp is False:
        os.system('rm %s.sam' % output_prefix)
        os.system('rm %s.bam' % output_prefix)


if __name__ == '__main__':

    reads2bam_parser = argparse.ArgumentParser(usage=reads2bam_usage)

    reads2bam_parser.add_argument('-p',               required=True,                                     help='output prefix')
    reads2bam_parser.add_argument('-ref',             required=True,                                     help='reference sequences')
    reads2bam_parser.add_argument('-index_ref',       required=False, action="store_true",               help='index reference')
    reads2bam_parser.add_argument('-r1',              required=False, default=None,                      help='paired reads r1')
    reads2bam_parser.add_argument('-r2',              required=False, default=None,                      help='paired reads r2')
    reads2bam_parser.add_argument('-u',               required=False, default=None,                      help='unpaired reads')
    reads2bam_parser.add_argument('-fastq',           required=False, action="store_true",               help='reads in fastq format')
    reads2bam_parser.add_argument('-t',               required=False, type=int, default=1,               help='number of threads, default: 1')
    reads2bam_parser.add_argument('-tmp',             required=False, action="store_true",               help='keep temporary files')

    args = vars(reads2bam_parser.parse_args())

    reads2bam(args)
