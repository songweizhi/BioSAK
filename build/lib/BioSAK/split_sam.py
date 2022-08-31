import os
import argparse
from time import sleep
import multiprocessing as mp


split_sam_usage = '''
==================== split_sam example command ====================

module load samtools/1.15
BioSAK split_sam -p single_ctg -i input.bam -r contig_1 -t 12
BioSAK split_sam -p multi_ctgs -i input.bam -r ctgs.txt -t 12

# Output files:
prefix.sam, prefix_sorted.bam and prefix_sorted.bam.bai

# ctgs.txt file format: one id per line, ">" excluded.

===================================================================
'''


def split_sam(args):

    sam_in        = args['i']
    input_ref_id  = args['r']
    output_prefix = args['p']
    thread_num    = args['t']

    # get input ref list
    input_ref_list = []
    if os.path.isfile(input_ref_id):
        for each_ref in open(input_ref_id):
            current_ref_id = each_ref.strip()
            input_ref_list.append(current_ref_id)
    else:
        input_ref_list = [input_ref_id]

    # get all header lines
    sam_all_header_lines_txt = '%s_all_header_lines.txt' % output_prefix
    get_sam_header_cmd = 'samtools view --header-only %s > %s' % (sam_in, sam_all_header_lines_txt)
    print(get_sam_header_cmd)
    os.system(get_sam_header_cmd)
    sleep(1)

    # define file name
    header_subset_txt    = '%s_header_lines.txt' % output_prefix
    sam_subset_no_header = '%s_no_header.sam'    % output_prefix
    subset_sam           = '%s.sam'              % output_prefix
    subset_bam           = '%s.bam'              % output_prefix
    subset_bam_sorted    = '%s_sorted.bam'       % output_prefix

    # get header line subset
    header_subset_txt_handle = open(header_subset_txt, 'w')
    for each_line in open(sam_all_header_lines_txt):
        if not each_line.startswith('@SQ'):
            header_subset_txt_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            ctg_id = each_line_split[1][3:]
            if ctg_id in input_ref_list:
                header_subset_txt_handle.write(each_line)
    header_subset_txt_handle.close()

    # extract reads aligned to provided references
    if not os.path.isfile(input_ref_id):

        get_aln_cmd = 'samtools view %s %s -o %s' % (sam_in, input_ref_id, sam_subset_no_header)
        print(get_aln_cmd)
        os.system(get_aln_cmd)
        sleep(1)

    else:
        split_sam_wd = '%s_split_sam_wd' % output_prefix
        os.mkdir(split_sam_wd)

        subset_sam_cmd_list = []
        for each_ref in open(input_ref_id):
            current_ref_id = each_ref.strip()
            get_aln_cmd = 'samtools view %s %s -o %s/%s.sam' % (sam_in, current_ref_id, split_sam_wd, current_ref_id)
            subset_sam_cmd_list.append(get_aln_cmd)

        # extract alignments
        pool = mp.Pool(processes=thread_num)
        pool.map(os.system, subset_sam_cmd_list)
        pool.close()
        pool.join()

        combine_extract_aln_cmd = 'cat %s/*.sam > %s' % (split_sam_wd, sam_subset_no_header)
        os.system(combine_extract_aln_cmd)

        # remove tmp dir
        os.system('rm -r %s' % split_sam_wd)

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


if __name__ == '__main__':

    split_sam_parser = argparse.ArgumentParser(usage=split_sam_usage)
    split_sam_parser.add_argument('-p', required=True, help='output prefix')
    split_sam_parser.add_argument('-i', required=True, help='input sam/bam file')
    split_sam_parser.add_argument('-r', required=True, help='reference id')
    split_sam_parser.add_argument('-t', required=False, type=int, default=1, help='number of threads')
    args = vars(split_sam_parser.parse_args())
    split_sam(args)
