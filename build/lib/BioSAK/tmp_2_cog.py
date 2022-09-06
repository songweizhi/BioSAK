
# file in
pwd_cog_20_def_tab      = '/Users/songweizhi/Desktop/cog-20.def.tab'
annotation_table        = '/Users/songweizhi/Desktop/eggnog_relab.tsv'

# file out
annotation_table_new    = '/Users/songweizhi/Desktop/eggnog_relab_with_description.tsv'


# get cog_id_to_description_dict (cognames2003-2014.tab)
cog_id_to_description_dict = {}
for cog_id_to_cate_des in open(pwd_cog_20_def_tab, encoding='windows-1252'):
    if not cog_id_to_cate_des.startswith('#'):
        cog_id_to_cate_des_split = cog_id_to_cate_des.strip().split('\t')
        cog_id = cog_id_to_cate_des_split[0]
        cog_cate = cog_id_to_cate_des_split[1]
        cog_des = cog_id_to_cate_des_split[2]
        cog_id_to_description_dict[cog_id] = cog_des


annotation_table_new_handle = open(annotation_table_new, 'w')
for each_line in open(annotation_table):
    if each_line.startswith('COG'):
        cog_id = each_line.strip().split('\t')[0]
        if '|' in cog_id:
            cog_id = cog_id.split('|')[0]
        cog_desc = cog_id_to_description_dict.get(cog_id, 'NA')
        annotation_table_new_handle.write('%s\t%s\n' % (each_line.strip(), cog_desc))
    else:
        annotation_table_new_handle.write('%s\t%s\n' % (each_line.strip(), 'NA'))
annotation_table_new_handle.close()

