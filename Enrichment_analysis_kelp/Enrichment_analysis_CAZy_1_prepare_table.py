import os
import glob

stats_GeneNumber_folder =   '/Users/songweizhi/Desktop/enrichment_analysis/1030_dbCAN_stats_GeneNumber'
file_out =                  '/Users/songweizhi/Desktop/enrichment_analysis/CAZy_enrichment_analysis.csv'

stats_GeneNumber_file_re = '%s/*_dbCAN_stats_GeneNumber.txt' % stats_GeneNumber_folder
stats_GeneNumber_file_list = [os.path.basename(file_name) for file_name in glob.glob(stats_GeneNumber_file_re)]


identified_cazy_list = []
cazy_dict_of_dict = {}
for cazy_stats_GeneNumber_file in stats_GeneNumber_file_list:

    MAG_id = cazy_stats_GeneNumber_file.split('_dbCAN_stats_GeneNumber.txt')[0]
    pwd_cog_stats_GeneNumber_file = '%s/%s' % (stats_GeneNumber_folder, cazy_stats_GeneNumber_file)

    current_MAG_cazy_stats_dict = {}
    current_MAG_cazy_total = 0
    for each in open(pwd_cog_stats_GeneNumber_file):

        if not each.startswith('Family	GeneNumber	Activities'):

            each_split = each.strip().split('\t')

            if each_split[0] not in identified_cazy_list:
                identified_cazy_list.append(each_split[0])

            current_MAG_cazy_stats_dict[each_split[0]] = int(each_split[1])
            current_MAG_cazy_total += int(each_split[1])

    current_MAG_cazy_stats_dict_norm = {}
    for cazy in current_MAG_cazy_stats_dict:
        current_MAG_cazy_stats_dict_norm[cazy] = float("{0:.3f}".format(current_MAG_cazy_stats_dict[cazy] * 100 / current_MAG_cazy_total))

    cazy_dict_of_dict[MAG_id] = current_MAG_cazy_stats_dict_norm


identified_cazy_list_sorted = sorted(identified_cazy_list)


file_out_handle = open(file_out, 'w')
file_out_handle.write('Source,MAG,%s\n' % (','.join(identified_cazy_list_sorted)))
for MAG in cazy_dict_of_dict:

    current_MAG_cazy_stats_list = []
    current_MAG_cog_stats = cazy_dict_of_dict[MAG]
    for cazy_id in identified_cazy_list_sorted:
        cazy_id_num = 0
        if cazy_id in current_MAG_cog_stats:
            cazy_id_num = current_MAG_cog_stats[cazy_id]
        current_MAG_cazy_stats_list.append(cazy_id_num)

    current_MAG_cazy_stats_list_str = [str(i) for i in current_MAG_cazy_stats_list]

    MAG_source = 'planktonic'
    if '_Refined_' in MAG:
        MAG_source = 'kelp-associated'

    file_out_handle.write('%s,%s,%s\n' % (MAG_source, MAG, ','.join(current_MAG_cazy_stats_list_str)))

file_out_handle.close()

