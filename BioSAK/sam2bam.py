import os
import argparse


sam2bam_usage = '''
============== sam2bam example command ==============

module load samtools
BioSAK sam2bam -sam input.sam -t 12

# This is a wrapper for the following three steps:
samtools view -bS input.sam -o input.bam
samtools sort input.bam -o input_sorted.bam
samtools index input_sorted.bam

=====================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def sam2bam(args):

    sam_in     = args['sam']
    thread_num = args['t']

    sam_path, sam_basename, sam_ext = sep_path_basename_ext(sam_in)

    pwd_bam             = '%s/%s.bam'                               % (sam_path, sam_basename)
    pwd_bam_sorted      = '%s/%s_sorted.bam'                        % (sam_path, sam_basename)
    cmd_samtools_view   = 'samtools view -bS %s -o %s --threads %s' % (sam_in, pwd_bam, thread_num)
    cmd_samtools_sort   = 'samtools sort %s -o %s --threads %s'     % (pwd_bam, pwd_bam_sorted, thread_num)
    cmd_samtools_index  = 'samtools index %s'                       % (pwd_bam_sorted)

    print('Converting sam to bam')
    os.system(cmd_samtools_view)
    print('Sorting bam file')
    os.system(cmd_samtools_sort)
    print('Indexing sorted bam file')
    os.system(cmd_samtools_index)
    print('Removing unsorted bam file')
    os.system('rm %s' % pwd_bam)
    print('All done!')

if __name__ == '__main__':

    sam2bam_parser = argparse.ArgumentParser(usage=sam2bam_usage)
    sam2bam_parser.add_argument('-sam', required=True,                       help='sam file')
    sam2bam_parser.add_argument('-t',   required=False, type=int, default=1, help='number of threads')
    args = vars(sam2bam_parser.parse_args())
    sam2bam(args)
