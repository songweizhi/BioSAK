
raw_file = '/Users/songweizhi/Desktop/DP_kegg.txt'

ko_id_to_num_dict = {}
for each in open(raw_file):
    each_split = each.strip().split('\t')
    if len(each_split) > 1:
        if each_split[1] not in ko_id_to_num_dict:
            ko_id_to_num_dict[each_split[1]] = 1
        else:
            ko_id_to_num_dict[each_split[1]] += 1

total_num = 0
for each_ko in ko_id_to_num_dict:
    total_num += ko_id_to_num_dict[each_ko]

print(total_num)


identified_ko = set()
for each in open('/Users/songweizhi/Desktop/DP_kegg_ko_stats_D_GeneNumber.txt'):
    each_split = each.strip().split('\t')
    #print(each_split)
    identified_ko.add(each_split[0])

for i in ko_id_to_num_dict:

    if i not in identified_ko:
        print(i)

'''
BioSAK COG2020 -m P -t 12 -db_dir /srv/scratch/z5039045/DB/COG2020 -i faa_files -x faa
'''