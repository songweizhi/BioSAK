
# file in/out
KEGG_DB_ko =                                        '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/ko00001.keg'
ko_level = 'D'
lifestyle = 'kelp'  # kelp tara
emriched_in_kelp_no_uncharacterized_protein =       '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_%s_no_uncharacterized_protein.txt'         % (ko_level, lifestyle)
emriched_in_kelp_no_uncharacterized_protein_ABCD =  '/Users/songweizhi/Desktop/Kelp_NM/enrichment_analysis/KEGG_%s_enrichment_test_emriched_in_%s_no_uncharacterized_protein_ABCD.txt'    % (ko_level, lifestyle)


# store ko functions in dict
As_description_dict = {}
Bs_description_dict = {}
Cs_description_dict = {}
Ds_description_dict = {}
D2ABCD_dict = {}
current_A = ''
current_B = ''
current_C = ''
for each_line in open(KEGG_DB_ko):
    if each_line[0] in ['A', 'B', 'C', 'D']:
        each_line_split = each_line.strip().split(' ')

        if each_line[0] == 'A':
            current_A_id = each_line_split[0]
            current_A_description = ' '.join(each_line_split[1:])
            current_A = current_A_id
            As_description_dict[current_A_id] = current_A_description

        elif each_line[0] == 'B':
            if len(each_line_split) > 1:
                current_B_id = each_line_split[2]
                current_B_description = ' '.join(each_line_split[3:])
                current_B = current_B_id
                Bs_description_dict[current_B_id] = current_B_description

        elif each_line[0] == 'C':
            current_C_id = each_line_split[4]
            current_C_description = ' '.join(each_line_split[5:])
            current_C = current_C_id
            Cs_description_dict[current_C_id] = current_C_description

        elif each_line[0] == 'D':
            current_D_id = each_line_split[6]
            current_D_description = ' '.join(each_line_split[7:])
            Ds_description_dict[current_D_id] = current_D_description
            ABCD_value = 'A_%s|B_%s|C_%s|D_%s' % (current_A, current_B, current_C, current_D_id)
            if current_D_id not in D2ABCD_dict:
                D2ABCD_dict[current_D_id] = [ABCD_value]
            elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                D2ABCD_dict[current_D_id].append(ABCD_value)


emriched_in_kelp_no_uncharacterized_protein_ABCD_handle = open(emriched_in_kelp_no_uncharacterized_protein_ABCD, 'w')
for each in open(emriched_in_kelp_no_uncharacterized_protein):

    if each.startswith('KO\tP_value'):
        emriched_in_kelp_no_uncharacterized_protein_ABCD_handle.write('%s_D\tKO_C\tFunction_C\tKO_B\tFunction_B\tKO_A\tFunction_A\n' % each.strip())
    else:
        ko_D = each.split('\t')[0]
        ko_D_to_ABCD = D2ABCD_dict[ko_D]
        current_ko_A_list = []
        current_ko_A_description_list = []
        current_ko_B_list = []
        current_ko_B_description_list = []
        current_ko_C_list = []
        current_ko_C_description_list = []
        for each_ABCD in ko_D_to_ABCD:
            each_ABCD_split = each_ABCD.split('|')
            if each_ABCD_split[0][2:] not in current_ko_A_list:
                current_ko_A_list.append(each_ABCD_split[0][2:])
                current_ko_A_description_list.append(As_description_dict[each_ABCD_split[0][2:]])
            if each_ABCD_split[1][2:] not in current_ko_B_list:
                current_ko_B_list.append(each_ABCD_split[1][2:])
                current_ko_B_description_list.append(Bs_description_dict[each_ABCD_split[1][2:]])
            if each_ABCD_split[2][2:] not in current_ko_C_list:
                current_ko_C_list.append(each_ABCD_split[2][2:])
                current_ko_C_description_list.append(Cs_description_dict[each_ABCD_split[2][2:]])

        emriched_in_kelp_no_uncharacterized_protein_ABCD_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each.strip(),
                                                                                                      '|'.join(current_ko_C_list),
                                                                                                      '|'.join(current_ko_C_description_list),
                                                                                                      '|'.join(current_ko_B_list),
                                                                                                      '|'.join(current_ko_B_description_list),
                                                                                                      '|'.join(current_ko_A_list),
                                                                                                      '|'.join(current_ko_A_description_list)))

emriched_in_kelp_no_uncharacterized_protein_ABCD_handle.close()

