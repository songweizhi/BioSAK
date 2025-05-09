import os
import argparse
import pandas as pd


rm_low_abd_otu_usage = '''
=================== rm_low_abd_otu example commands ===================

# ignore OTUs with relative abundance lower than 0.1%
BioSAK rm_low_abd_otu -i otu_table.txt -c 0.001 -o otu_table_0.001.txt

# ignore OTUs with relative abundance lower than 1%
BioSAK rm_low_abd_otu -i otu_table.txt -c 0.01 -o otu_table_0.01.txt

=======================================================================
'''


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):
    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def rm_low_abd_otu(args):

    otu_table_in    = args['i']
    abd_cutoff      = args['c']
    otu_table_out   = args['o']

    # define tmp file name
    tmp1_t          = '%s.tmp1' % otu_table_out
    tmp2_t_filtered = '%s.tmp2' % otu_table_out

    transpose_csv(otu_table_in, tmp1_t, '\t', 0, 0)

    otu_table_txt_t_filtered_handle = open(tmp2_t_filtered, 'w')
    line_index = 0
    for each_line in open(tmp1_t):
        if line_index == 0:
            otu_table_txt_t_filtered_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            sample_id = each_line_split[0]
            count_list = [int(i) for i in each_line_split[1:]]
            count_sum = sum(count_list)
            count_list_filtered = []
            for each_count in count_list:
                if (each_count/count_sum) >= abd_cutoff:
                    count_list_filtered.append(str(each_count))
                else:
                    count_list_filtered.append('0')
            otu_table_txt_t_filtered_handle.write('%s\t%s\n' % (sample_id, '\t'.join(count_list_filtered)))
        line_index += 1
    otu_table_txt_t_filtered_handle.close()

    transpose_csv(tmp2_t_filtered, otu_table_out, '\t', 0, 0)

    # remove tmp files
    os.remove(tmp1_t)
    os.remove(tmp2_t_filtered)
    print('Done')


if __name__ == '__main__':

    rm_low_abd_otu_parser = argparse.ArgumentParser()
    rm_low_abd_otu_parser.add_argument('-i', required=True,                              help='input otu count table')
    rm_low_abd_otu_parser.add_argument('-c', required=False, default=0.001,type=float,   help='relative abundance cutoff, default is 0.001')
    rm_low_abd_otu_parser.add_argument('-o', required=True,                              help='output otu count table')
    args = vars(rm_low_abd_otu_parser.parse_args())
    rm_low_abd_otu(args)
