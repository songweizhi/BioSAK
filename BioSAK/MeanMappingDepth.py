import argparse
import numpy as np


MeanMappingDepth_usage = '''
========== MeanMappingDepth example commands ==========

BioSAK MeanMappingDepth -depth ctgs.bam.depth
BioSAK MeanMappingDepth -depth ctgs.bam.depth -T

# get depth file with samtools
samtools depth ctgs.bam > ctgs.bam.depth

=======================================================
'''


def MeanMappingDepth(args):

    sam_depth  = args['depth']
    total_mean = args['T']

    # get length and total depth
    seq_to_len_and_depth_dict = {}
    for each_line in open(sam_depth):
        each_line_split = each_line.strip().split('\t')
        seq_id = each_line_split[0]
        bp_depth = int(each_line_split[2])
        if seq_id not in seq_to_len_and_depth_dict:
            seq_to_len_and_depth_dict[seq_id] = [1, bp_depth]
        else:
            seq_to_len_and_depth_dict[seq_id][0] += 1
            seq_to_len_and_depth_dict[seq_id][1] += bp_depth

    # get mean depth
    all_seq_len   = 0
    all_seq_depth = 0
    all_seq_mean_depth_list = []
    if total_mean is False:
        print('Seq\tLength(bp)\tMeanDepth')
    for each_seq in seq_to_len_and_depth_dict:
        current_seq_len   = seq_to_len_and_depth_dict[each_seq][0]
        current_seq_depth = seq_to_len_and_depth_dict[each_seq][1]
        all_seq_len += current_seq_len
        all_seq_depth += current_seq_depth
        seq_mean_depth = float("{0:.2f}".format(current_seq_depth/current_seq_len))
        if total_mean is False:
            print('%s\t%s\t%s' % (each_seq, seq_to_len_and_depth_dict[each_seq][0], seq_mean_depth))

        all_seq_mean_depth_list.append(seq_mean_depth)

    if total_mean is True:
        all_seq_mean_depth = float("{0:.2f}".format(all_seq_depth/all_seq_len))

        # report
        print('Sequence num  = %s'      % len(seq_to_len_and_depth_dict))
        print('Total length  = %s bp'   % all_seq_len)
        print('Highest depth = %s'      % max(all_seq_mean_depth_list))
        print('Lowest depth  = %s'      % min(all_seq_mean_depth_list))
        print('Median depth  = %s'      % np.median(all_seq_mean_depth_list))
        print('Mean depth    = %s'      % all_seq_mean_depth)


if __name__ == '__main__':

    MeanMappingDepth_parser = argparse.ArgumentParser()

    MeanMappingDepth_parser.add_argument('-depth', required=True,                       help='input depth file from "samtools depth" ')
    MeanMappingDepth_parser.add_argument('-T',     required=False, action="store_true", help='get overall stats')

    args = vars(MeanMappingDepth_parser.parse_args())

    MeanMappingDepth(args)
