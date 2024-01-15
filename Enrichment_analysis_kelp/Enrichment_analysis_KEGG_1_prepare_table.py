import os
import glob


ko_level = 'B'  # C D
stats_GeneNumber_folder =   '/Users/songweizhi/Desktop/enrichment_analysis/1030_ko_stats_%s_GeneNumber'     % ko_level
file_out =                  '/Users/songweizhi/Desktop/enrichment_analysis/KEGG_%s_enrichment_analysis.csv' % ko_level


stats_GeneNumber_file_re = '%s/*_ko_stats_%s_GeneNumber.txt' % (stats_GeneNumber_folder, ko_level)
stats_GeneNumber_file_list = [os.path.basename(file_name) for file_name in glob.glob(stats_GeneNumber_file_re)]

identified_ko_list = []
ko_dict_of_dict = {}
for ko_stats_GeneNumber_file in stats_GeneNumber_file_list:

    MAG_id = ko_stats_GeneNumber_file.split('_ko_stats_%s_GeneNumber.txt' % ko_level)[0]
    pwd_ko_stats_GeneNumber_file = '%s/%s' % (stats_GeneNumber_folder, ko_stats_GeneNumber_file)

    current_MAG_ko_stats_dict = {}
    current_MAG_ko_total = 0
    for each in open(pwd_ko_stats_GeneNumber_file):

        if not each.startswith('KO	GeneNumber	Description'):

            each_split = each.strip().split('\t')

            if each_split[0] not in identified_ko_list:
                identified_ko_list.append(each_split[0])

            current_MAG_ko_stats_dict[each_split[0]] = int(each_split[1])
            current_MAG_ko_total += int(each_split[1])

    current_MAG_ko_stats_dict_norm = {}
    for cog in current_MAG_ko_stats_dict:
        current_MAG_ko_stats_dict_norm[cog] = float("{0:.3f}".format(current_MAG_ko_stats_dict[cog] * 100 / current_MAG_ko_total))

    ko_dict_of_dict[MAG_id] = current_MAG_ko_stats_dict_norm


identified_ko_list_sorted = sorted(identified_ko_list)


file_out_handle = open(file_out, 'w')
file_out_handle.write('Source,MAG,%s\n' % (','.join(identified_ko_list_sorted)))
for MAG in ko_dict_of_dict:

    current_MAG_ko_stats_list = []
    current_MAG_ko_stats = ko_dict_of_dict[MAG]
    for ko_id in identified_ko_list_sorted:
        ko_id_num = 0
        if ko_id in current_MAG_ko_stats:
            ko_id_num = current_MAG_ko_stats[ko_id]
        current_MAG_ko_stats_list.append(ko_id_num)

    current_MAG_ko_stats_list_str = [str(i) for i in current_MAG_ko_stats_list]

    MAG_source = 'planktonic'
    if '_Refined_' in MAG:
        MAG_source = 'kelp-associated'

    file_out_handle.write('%s,%s,%s\n' % (MAG_source, MAG, ','.join(current_MAG_ko_stats_list_str)))

file_out_handle.close()

