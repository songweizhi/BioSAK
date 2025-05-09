import argparse
import pandas as pd


rm_low_depth_sample_usage = '''
===================== rm_low_depth_sample example commands =====================

BioSAK rm_low_depth_sample -i otu_table.txt -c 20000 -o otu_table_mim20000.txt

================================================================================
'''


def rm_low_depth_sample(args):

    otu_table_in            = args['i']
    abd_cutoff              = args['c']
    otu_table_out           = args['o']
    otu_table_df            = pd.read_csv(otu_table_in, sep='\t', header=0, index_col=0)
    column_sums             = otu_table_df.sum()
    column_sum_dict         = column_sums.to_dict()
    otu_table_df_filtered   = otu_table_df.loc[:, column_sums >= abd_cutoff]
    otu_table_df_filtered.to_csv(otu_table_out, sep='\t')

    # report
    print('The following samples were removed from the OTU table:')
    for sample in sorted(column_sum_dict.keys()):
        seq_count = column_sum_dict[sample]
        if seq_count < abd_cutoff:
            print('%s\t%s' % (sample, seq_count))
    print('Done!')


if __name__ == '__main__':

    rm_low_depth_sample_parser = argparse.ArgumentParser()
    rm_low_depth_sample_parser.add_argument('-i', required=True,                            help='input otu count table')
    rm_low_depth_sample_parser.add_argument('-c', required=False, default=20000,type=int,   help='minimal number of sequences, default is 20000')
    rm_low_depth_sample_parser.add_argument('-o', required=True,                            help='output otu count table')
    args = vars(rm_low_depth_sample_parser.parse_args())
    rm_low_depth_sample(args)
