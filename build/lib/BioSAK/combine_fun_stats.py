import argparse


combine_fun_stats_usage = '''
========================================= combine_fun_stats example commands =========================================

BioSAK combine_fun_stats -f1 sponge_ko_c_stats.txt -f2 seawater_ko_c_stats.txt -l1 Sponge -l2 Seawater -o summary.txt

======================================================================================================================
'''

def combine_fun_stats(args):

    stats_1_txt     = args['f1']
    stats_1_label   = args['l1']
    stats_2_txt     = args['f2']
    stats_2_label   = args['l2']
    op_txt          = args['o']

    overall_fun_set = set()
    fun_desc_dict = dict()
    stats_1_dict = dict()
    stats_1_sum = 0
    for each_line in open(stats_1_txt):
        each_line_split = each_line.strip().split('\t')
        fun_id = each_line_split[0]
        fun_num = int(each_line_split[1])
        fun_desc = each_line_split[2]
        stats_1_dict[fun_id] = fun_num
        fun_desc_dict[fun_id] = fun_desc
        stats_1_sum += fun_num
        overall_fun_set.add(fun_id)

    stats_2_dict = dict()
    stats_2_sum = 0
    for each_line in open(stats_2_txt):
        each_line_split = each_line.strip().split('\t')
        fun_id = each_line_split[0]
        fun_num = int(each_line_split[1])
        fun_desc = each_line_split[2]
        stats_2_dict[fun_id] = fun_num
        fun_desc_dict[fun_id] = fun_desc
        stats_2_sum += fun_num
        overall_fun_set.add(fun_id)

    stats_1_dict_pct = dict()
    for i in stats_1_dict:
        i_pct = stats_1_dict[i]*100/stats_1_sum
        i_pct = float("{0:.2f}".format(i_pct))
        stats_1_dict_pct[i] = i_pct

    stats_2_dict_pct = dict()
    for i in stats_2_dict:
        i_pct = stats_2_dict[i]*100/stats_2_sum
        i_pct = float("{0:.2f}".format(i_pct))
        stats_2_dict_pct[i] = i_pct

    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('\t%s\t%s\tRatio\n' % (stats_1_label, stats_2_label))
    for each_fun in sorted(list(overall_fun_set)):
        fun_desc = fun_desc_dict[each_fun]
        fun_pct_1 = stats_1_dict_pct.get(each_fun, 0)
        fun_pct_2 = stats_2_dict_pct.get(each_fun, 0)

        ratio = 9999
        if fun_pct_2 != 0:
            ratio = fun_pct_1/fun_pct_2
        ratio = float("{0:.2f}".format(ratio))

        op_txt_handle.write('%s__%s\t%s\t%s\t%s\n' % (each_fun, fun_desc, fun_pct_1, fun_pct_2, ratio))
    op_txt_handle.close()


if __name__ == '__main__':

    combine_fun_stats_parser = argparse.ArgumentParser(usage=combine_fun_stats_usage)
    combine_fun_stats_parser.add_argument('-f1',    required=True,  help='fun stats 1')
    combine_fun_stats_parser.add_argument('-f2',    required=True,  help='fun stats 2')
    combine_fun_stats_parser.add_argument('-l1',    required=True,  help='fun stats 1 label')
    combine_fun_stats_parser.add_argument('-l2',    required=True,  help='fun stats 2 label')
    combine_fun_stats_parser.add_argument('-o',     required=True,  help='output directory')
    args = vars(combine_fun_stats_parser.parse_args())
    combine_fun_stats(args)
