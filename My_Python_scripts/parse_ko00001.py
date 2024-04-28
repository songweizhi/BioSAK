
KEGG_DB = 'ko00001.keg'

As_description_dict = {}
Bs_description_dict = {}
Cs_description_dict = {}
Ds_description_dict = {}
D2ABCD_dict = {}
current_A = ''
current_B = ''
current_C = ''
for each_line in open(KEGG_DB):
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

