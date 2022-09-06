
# file in
KEGG_DB_ko          = '/Users/songweizhi/Desktop/ko00001.keg'
annotation_table    = '/Users/songweizhi/Desktop/ko_relab_renamed.tsv'

# file out
annotation_table_new    = '/Users/songweizhi/Desktop/ko_relab_renamed_with_description.tsv'


Ds_description_dict = {}
for each_line in open(KEGG_DB_ko):
    if each_line[0] == 'D':
        each_line_split = each_line.strip().split(' ')
        current_D_id = each_line_split[6]
        current_D_description = ' '.join(each_line_split[7:])
        Ds_description_dict[current_D_id] = current_D_description

annotation_table_new_handle = open(annotation_table_new, 'w')
for each_line in open(annotation_table):
    if each_line.startswith('K'):
        ko_id = each_line.strip().split(':')[0]
        ko_desc = Ds_description_dict.get(ko_id, 'NA')
        annotation_table_new_handle.write('%s\t%s\n' % (each_line.strip(), ko_desc))
    else:
        annotation_table_new_handle.write('%s\t%s\n' % (each_line.strip(), 'NA'))
annotation_table_new_handle.close()
