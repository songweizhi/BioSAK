
def get_ko2description_dict(ko00001_keg):

    As_description_dict = {}
    Bs_description_dict = {}
    Cs_description_dict = {}
    Ds_description_dict = {}
    D2ABCD_dict = {}
    ko2level_dict = {}

    current_A = ''
    current_B = ''
    current_C = ''
    for each_line in open(ko00001_keg):
        if each_line[0] in ['A', 'B', 'C', 'D']:
            each_line_split = each_line.strip().split(' ')

            if each_line[0] == 'A':
                current_A_id = each_line_split[0]
                current_A_description = ' '.join(each_line_split[1:-1])
                ko2level_dict[current_A_id] = 'A'
                current_A = current_A_id
                As_description_dict[current_A_id] = current_A_description

            elif each_line[0] == 'B':
                if len(each_line_split) > 1:
                    current_B_id = each_line_split[2]
                    current_B_description = ' '.join(each_line_split[3:-1])
                    ko2level_dict[current_B_id] = 'B'
                    current_B = current_B_id
                    Bs_description_dict[current_B_id] = current_B_description

            elif each_line[0] == 'C':
                current_C_id = each_line_split[4]
                current_C_description = ' '.join(each_line_split[5:-1])
                ko2level_dict[current_C_id] = 'C'
                current_C = current_C_id
                Cs_description_dict[current_C_id] = current_C_description

            elif each_line[0] == 'D':
                current_D_id = each_line_split[6]
                current_D_description = ' '.join(each_line_split[7:])
                ko2level_dict[current_D_id] = 'D'
                Ds_description_dict[current_D_id] = current_D_description
                ABCD_value = '%s|%s|%s|%s' % (current_A, current_B, current_C, current_D_id)
                if current_D_id not in D2ABCD_dict:
                    D2ABCD_dict[current_D_id] = [ABCD_value]
                elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                    D2ABCD_dict[current_D_id].append(ABCD_value)

    return As_description_dict, Bs_description_dict, Cs_description_dict, Ds_description_dict, D2ABCD_dict, ko2level_dict


def get_ec_of_interested_ko(D2ABCD_dict, KO_description_D_dict, ko_level, ko_id):

    # get ec list
    interested_ko_ec_list = set()
    for each_ko_D in D2ABCD_dict:

        each_ko_D_description = KO_description_D_dict[each_ko_D]
        if '[EC:' in each_ko_D_description:
            ec_str = each_ko_D_description.strip().split('[EC:')[1][:-1]

            # get ec list
            ec_list = [ec_str]
            if ' ' in ec_str:
                ec_list = ec_str.split(' ')

            for each_group in D2ABCD_dict[each_ko_D]:
                each_ko_D_split = each_group.split('|')
                if (ko_level == 'A') and (each_ko_D_split[0] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'B') and (each_ko_D_split[1] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'C') and (each_ko_D_split[2] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

                elif (ko_level == 'D') and (each_ko_D_split[3] == ko_id):
                    for ec in ec_list:
                        interested_ko_ec_list.add(ec)

    return interested_ko_ec_list


ko00001_keg =      '/Users/songweizhi/DB/KEGG_2016-09-26/ko00001.keg'
interested_ko_id = '00010'


# read in KEGG db file
KO_description_A_dict, KO_description_B_dict, KO_description_C_dict, KO_description_D_dict, D2ABCD_dict, ko2level_dict = get_ko2description_dict(ko00001_keg)

# get ec list from interested KO category
interested_ec_list = get_ec_of_interested_ko(D2ABCD_dict, KO_description_D_dict, ko2level_dict[interested_ko_id], interested_ko_id)


print(interested_ec_list)
print(len(interested_ec_list))

