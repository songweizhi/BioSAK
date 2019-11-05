import os

for each_file in open('/srv/scratch/z5039045/MetaCHIP_demo/TT_90MGs/GoodBins_0.5_0.05_MetaCHIP_wd/faa_files.txt'):
    each_file_no_ext = each_file.strip().split('.')[0]

    cmd = 'cp %s_COG_results/func_stats.txt 0func_stats_files/%s_func_stats.txt' % (each_file_no_ext, each_file_no_ext)

    os.system(cmd)


