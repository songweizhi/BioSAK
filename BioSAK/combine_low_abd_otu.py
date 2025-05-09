import os
import argparse
import pandas as pd


combine_low_abd_otu_usage = '''
================== combine_low_abd_otu example commands ==================

# keep OTUs with relative abundance higher than 50% and combine the rests
BioSAK combine_low_abd_otu -i otu_table.txt -c 0.5 -o otu_table_combine_min_50.txt

==========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):
    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def combine_low_abd_otu(args):

    otu_table_in          = args['i']
    abd_cutoff            = args['c']
    otu_table_out         = args['o']
    interested_sample_txt = args['s']

    interested_sample_set = set()
    if interested_sample_txt is not None:
        if os.path.isfile(interested_sample_txt):
            for sample in open(interested_sample_txt):
                interested_sample_set.add(sample.strip().split()[0])
        else:
            print('%s not found, program exited!' % interested_sample_txt)
            exit()

    f_name, f_path, f_base, f_ext = sep_path_basename_ext(otu_table_out)

    # define tmp file name
    tmp1_t               = '%s.tmp1'        % otu_table_out
    tmp2_t_filtered      = '%s.tmp2'        % otu_table_out
    tmp2_t_filtered_norm = '%s.tmp2.norm'   % otu_table_out
    otu_table_out_norm   = '%s/%s_pct.%s'   % (f_path, f_base, f_ext)

    transpose_csv(otu_table_in, tmp1_t, '\t', 0, 0)

    # get df sample set
    df_sample_set = set()
    line_index = 0
    for each_line in open(tmp1_t):
        each_line_split = each_line.strip().split('\t')
        if line_index > 0:
            sample_id = each_line_split[0]
            df_sample_set.add(sample_id)

    if len(interested_sample_set) == 0:
        interested_sample_set = df_sample_set

    # get all OTUs that have abundance higher than specified cutoff in at least one sample
    otus_to_keep = set()
    line_index = 0
    otu_id_list = []
    for each_line in open(tmp1_t):
        each_line_split = each_line.strip().split('\t')
        if line_index == 0:
            otu_id_list = each_line_split
        else:
            sample_id = each_line_split[0]
            if sample_id in interested_sample_set:
                count_list = [int(i) for i in each_line_split[1:]]
                count_sum = sum(count_list)
                for (otu_id, otu_count) in zip(otu_id_list, count_list):
                    if (otu_count / count_sum) >= abd_cutoff:
                        otus_to_keep.add(otu_id)
        line_index += 1

    # filter
    otu_count_dict = dict()
    line_index = 0
    otu_id_list = []
    for each_line in open(tmp1_t):
        each_line_split = each_line.strip().split('\t')
        if line_index == 0:
            otu_id_list = each_line_split
        else:
            sample_id = each_line_split[0]
            if sample_id in interested_sample_set:
                otu_count_list = [int(i) for i in each_line_split[1:]]
                current_sample_otu_count_dict = dict()
                for (otu_id, otu_count) in zip(otu_id_list, otu_count_list):
                    if otu_id in otus_to_keep:
                        current_sample_otu_count_dict[otu_id] = otu_count
                    else:
                        if 'others' not in current_sample_otu_count_dict:
                            current_sample_otu_count_dict['others'] = 0
                        current_sample_otu_count_dict['others'] += otu_count
                otu_count_dict[sample_id] = current_sample_otu_count_dict
        line_index += 1

    otus_to_keep.add('others')
    otus_to_keep_list_sorted = sorted(list(otus_to_keep))

    # write out
    otu_table_txt_t_filtered_handle = open(tmp2_t_filtered, 'w')
    otu_table_txt_t_filtered_handle.write('\t%s\n' % '\t'.join(otus_to_keep_list_sorted))
    for each_sample in otu_count_dict:
        combined_otu_count_dict = otu_count_dict[each_sample]
        count_list= [each_sample]
        for eahc_otu in otus_to_keep_list_sorted:
            count_list.append(combined_otu_count_dict[eahc_otu])
        count_list = [str(i) for i in count_list]
        otu_table_txt_t_filtered_handle.write('\t'.join(count_list) + '\n')
    otu_table_txt_t_filtered_handle.close()

    tmp2_t_filtered_norm_handle = open(tmp2_t_filtered_norm, 'w')
    line_index = 0
    for each_line in open(tmp2_t_filtered):
        each_line_split = each_line.strip().split('\t')
        if line_index == 0:
            tmp2_t_filtered_norm_handle.write(each_line)
        else:
            sample_id = each_line_split[0]
            value_list = [int(i) for i in each_line_split[1:]]
            otu_sum = sum(value_list)
            value_list_norm = [float("{0:.4f}".format(i*100/otu_sum)) for i in value_list]
            tmp2_t_filtered_norm_handle.write('%s\t%s\n' % (sample_id, '\t'.join([str(i) for i in value_list_norm])))
        line_index += 1
    tmp2_t_filtered_norm_handle.close()

    transpose_csv(tmp2_t_filtered, otu_table_out, '\t', 0, 0)
    transpose_csv(tmp2_t_filtered_norm, otu_table_out_norm, '\t', 0, 0)


    # remove tmp files
    os.remove(tmp1_t)
    os.remove(tmp2_t_filtered)
    os.remove(tmp2_t_filtered_norm)
    print('Results exported to:\n%s\n%s' % (otu_table_out, otu_table_out_norm))
    print('Done!')


if __name__ == '__main__':

    combine_low_abd_otu_parser = argparse.ArgumentParser()
    combine_low_abd_otu_parser.add_argument('-i', required=True,                            help='input otu count table')
    combine_low_abd_otu_parser.add_argument('-c', required=False, default=0.5,type=float,   help='relative abundance cutoff, default is 0.5')
    combine_low_abd_otu_parser.add_argument('-s', required=False, default=None,             help='interested sample txt')
    combine_low_abd_otu_parser.add_argument('-o', required=True,                            help='output otu table')
    args = vars(combine_low_abd_otu_parser.parse_args())
    combine_low_abd_otu(args)


'''

python3 /Users/songweizhi/PycharmProjects/BioSAK/BioSAK/combine_low_abd_otu.py -i /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250325/s07_AllSamples_unoise_otu_table_noEU_mim20000.txt -o /Users/songweizhi/Desktop/SMP/02_Usearch_BLCA_GTDB_20250325/s07_AllSamples_unoise_otu_table_noEU_mim20000_combined_low_abundant_OTUs_0.5.txt -c 0.1

'''
