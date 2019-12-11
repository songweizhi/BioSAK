

for each in open('/Users/songweizhi/Desktop/good_bins_id.txt'):

    cp_cmd = 'cp MetaBAT_bins_Test_c17_faa_files/%s_COG_results/func_stats.txt func_stats_txt_files/%s_func_stats.txt' % (each.strip(), each.strip())
    print(cp_cmd)

