import os
import glob

# what does this script do:
# 1. transfer absolute number into percentage
# 2. prepare table

cog_stats_GeneNumber_folder = '/Users/songweizhi/Desktop/1030_cog_stats_GeneNumber'
file_out = '/Users/songweizhi/Desktop/COG_enrichment_analysis.csv'


cog_stats_GeneNumber_file_re = '%s/*_cog_stats_GeneNumber.txt' % cog_stats_GeneNumber_folder
cog_stats_GeneNumber_file_list = [os.path.basename(file_name) for file_name in glob.glob(cog_stats_GeneNumber_file_re)]


identified_cog_list = []
cog_dict_of_dict = {}
for cog_stats_GeneNumber_file in cog_stats_GeneNumber_file_list:

    MAG_id = cog_stats_GeneNumber_file.split('_cog_stats_GeneNumber')[0]
    pwd_cog_stats_GeneNumber_file = '%s/%s' % (cog_stats_GeneNumber_folder, cog_stats_GeneNumber_file)

    current_MAG_cog_stats_dict = {}
    current_MAG_cog_total = 0
    for each in open(pwd_cog_stats_GeneNumber_file):

        if not each.startswith('COG	GeneNumber	Description'):

            each_split = each.strip().split('\t')

            if each_split[0] not in identified_cog_list:
                identified_cog_list.append(each_split[0])

            current_MAG_cog_stats_dict[each_split[0]] = int(each_split[1])
            current_MAG_cog_total += int(each_split[1])

    current_MAG_cog_stats_dict_norm = {}
    for cog in current_MAG_cog_stats_dict:
        current_MAG_cog_stats_dict_norm[cog] = float("{0:.3f}".format(current_MAG_cog_stats_dict[cog]*100/current_MAG_cog_total))

    cog_dict_of_dict[MAG_id] = current_MAG_cog_stats_dict_norm
    #cog_dict_of_dict[MAG_id] = current_MAG_cog_stats_dict


identified_cog_list_sorted = sorted(identified_cog_list)


file_out_handle = open(file_out, 'w')
file_out_handle.write('Source,MAG,%s\n' % (','.join(identified_cog_list_sorted)))
for MAG in cog_dict_of_dict:

    current_MAG_cog_stats_list = []
    current_MAG_cog_stats = cog_dict_of_dict[MAG]
    for cog_id in identified_cog_list_sorted:
        cog_id_num = 0
        if cog_id in current_MAG_cog_stats:
            cog_id_num = current_MAG_cog_stats[cog_id]
        current_MAG_cog_stats_list.append(cog_id_num)

    current_MAG_cog_stats_list_str = [str(i) for i in current_MAG_cog_stats_list]

    MAG_source = 'planktonic'
    if '_Refined_' in MAG:
        MAG_source = 'kelp-associated'

    file_out_handle.write('%s,%s,%s\n' % (MAG_source, MAG, ','.join(current_MAG_cog_stats_list_str)))

file_out_handle.close()

