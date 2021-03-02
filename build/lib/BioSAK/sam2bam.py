import os
import argparse
from BioSAK.global_functions import sep_path_basename_ext


sam2bam_usage = '''
========== sam2bam example command ==========

module load samtools/1.10
BioSAK sam2bam -sam input.sam

=============================================
'''


def sam2bam(args):

    sam_in   = args['sam']

    sam_path, sam_basename, sam_ext = sep_path_basename_ext(sam_in)

    pwd_bam         = '%s/%s.bam'           % (sam_path, sam_basename)
    pwd_bam_sorted  = '%s/%s_sorted.bam'    % (sam_path, sam_basename)


    cmd_samtools_view   = 'samtools view -bS %s -o %s'  % (sam_in, pwd_bam)
    cmd_samtools_sort   = 'samtools sort %s -o %s'      % (pwd_bam, pwd_bam_sorted)
    cmd_samtools_index  = 'samtools index %s'           % pwd_bam_sorted

    os.system(cmd_samtools_view)
    os.system(cmd_samtools_sort)
    os.system(cmd_samtools_index)

    os.system('rm %s' % pwd_bam)


if __name__ == '__main__':

    sam2bam_parser = argparse.ArgumentParser(usage=sam2bam_usage)
    sam2bam_parser.add_argument('-sam', required=True, help='sam file')
    args = vars(sam2bam_parser.parse_args())
    sam2bam(args)
