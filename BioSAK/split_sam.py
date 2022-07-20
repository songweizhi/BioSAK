import os
from time import sleep
import argparse


split_sam_usage = '''
================== split_sam example command ==================

module load samtools/1.15
BioSAK split_sam -p contig_1 -i input.sam -r contig_1 -t 12
BioSAK split_sam -p contig_1 -i input.bam -r contig_1 -t 12

# Output files:
contig_1.sam
contig_1_sorted.bam
contig_1_sorted.bam.bai

===============================================================
'''


def split_sam(args):

    sam_in        = args['i']
    ref_id        = args['r']
    output_prefix = args['p']
    thread_num    = args['t']

    sam_all_header_lines_txt = '%s_all_header_lines.txt' % output_prefix
    get_sam_header_cmd = 'samtools view --header-only %s > %s' % (sam_in, sam_all_header_lines_txt)
    print(get_sam_header_cmd)
    os.system(get_sam_header_cmd)
    sleep(1)

    if not os.path.isfile(ref_id):

        sam_subset_no_header = '%s_no_header.sam'    % output_prefix
        header_subset_txt    = '%s_header_lines.txt' % output_prefix
        subset_sam           = '%s.sam'              % output_prefix
        subset_bam           = '%s.bam'              % output_prefix
        subset_bam_sorted    = '%s_sorted.bam'       % output_prefix

        get_aln_cmd = 'samtools view %s %s -o %s' % (sam_in, ref_id, sam_subset_no_header)
        print(get_aln_cmd)
        os.system(get_aln_cmd)
        sleep(1)

        header_subset_txt_handle = open(header_subset_txt, 'w')
        for each_line in open(sam_all_header_lines_txt):
            if not each_line.startswith('@SQ'):
                header_subset_txt_handle.write(each_line)
            else:
                each_line_split = each_line.strip().split('\t')
                ctg_id = each_line_split[1][3:]
                if ctg_id == ref_id:
                    header_subset_txt_handle.write(each_line)
        header_subset_txt_handle.close()

        # combine header and aligned reads
        os.system('cat %s %s > %s' % (header_subset_txt, sam_subset_no_header, subset_sam))
        sleep(1)

        # remove tmp files
        os.system('rm %s' % sam_all_header_lines_txt)
        os.system('rm %s' % sam_subset_no_header)
        os.system('rm %s' % header_subset_txt)

        # sam to bam
        cmd_samtools_view  = 'samtools view -bS %s -o %s --threads %s' % (subset_sam, subset_bam, thread_num)
        print(cmd_samtools_view)
        os.system(cmd_samtools_view)
        sleep(1)

        # sort bam
        cmd_samtools_sort  = 'samtools sort %s -o %s --threads %s' % (subset_bam, subset_bam_sorted, thread_num)
        print(cmd_samtools_sort)
        os.system(cmd_samtools_sort)
        sleep(1)

        # index bam
        cmd_samtools_index = 'samtools index %s' % subset_bam_sorted
        print(cmd_samtools_index)
        os.system(cmd_samtools_index)
        sleep(1)

        print('Removing unsorted bam file')
        os.system('rm %s' % subset_bam)
        print('Done!')

    else:
        pass  # to be added (for multiple reference sequences)


if __name__ == '__main__':

    split_sam_parser = argparse.ArgumentParser(usage=split_sam_usage)
    split_sam_parser.add_argument('-p', required=True, help='output prefix')
    split_sam_parser.add_argument('-i', required=True, help='input sam/bam file')
    split_sam_parser.add_argument('-r', required=True, help='reference id')
    split_sam_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads')
    args = vars(split_sam_parser.parse_args())
    split_sam(args)
