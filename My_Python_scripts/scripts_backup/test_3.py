
pwd_cog2003_2014 = '/Users/songweizhi/Desktop/cog2003-2014.csv'

protein_id_to_cog_id_dict = {}
for protein_to_cog in open(pwd_cog2003_2014):
    protein_to_cog_split = protein_to_cog.strip().split(',')
    protein_id = protein_to_cog_split[2]
    cog_id = protein_to_cog_split[6]
    if protein_id not in protein_id_to_cog_id_dict:
        protein_id_to_cog_id_dict[protein_id] = {cog_id}
    else:
        protein_id_to_cog_id_dict[protein_id].add(cog_id)




for each_protein in protein_id_to_cog_id_dict:
    if len(protein_id_to_cog_id_dict[each_protein]) > 1:
        print('%s\t%s' % (each_protein, protein_id_to_cog_id_dict[each_protein]))


