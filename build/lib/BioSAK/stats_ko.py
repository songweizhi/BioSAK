import os
import argparse


stats_ko_usage = '''
=================================== stats_ko example commands ===================================

BioSAK stats_ko -ko ko_d.txt -db ko00001.keg -p Demo

# Required DB files:
ko00001.keg: https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext&filedir=
=================================================================================================
'''


def ko_stats_to_txt(ko_stats_dict, ko_desc_dict, op_txt, op_stats_txt):

    op_txt_handle = open(op_txt, 'w')
    op_stats_txt_handle = open(op_stats_txt, 'w')
    for ko_high in sorted(list(ko_stats_dict.keys())):
        ko_d_set = ko_stats_dict[ko_high]
        ko_high_desc = ko_desc_dict[ko_high]
        op_stats_txt_handle.write('%s\t%s\t%s\n' % (ko_high, len(ko_d_set), ko_high_desc))
        for ko_d in sorted(list(ko_d_set)):
            ko_d_desc = ko_desc_dict[ko_d]
            op_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (ko_high, len(ko_d_set), ko_high_desc.strip(), ko_d, ko_d_desc.strip()))
    op_txt_handle.close()
    op_stats_txt_handle.close()


def stats_ko(args):

    ko_txt              = args['ko']
    op_dir              = args['o']
    db_file             = args['db']
    force_create_op_dir = args['f']

    ####################################################################################################################

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    op_a_txt        = '%s/ko_A.txt'       % op_dir
    op_a_stats_txt  = '%s/ko_A_stats.txt' % op_dir
    op_b_txt        = '%s/ko_B.txt'       % op_dir
    op_b_stats_txt  = '%s/ko_B_stats.txt' % op_dir
    op_c_txt        = '%s/ko_C.txt'       % op_dir
    op_c_stats_txt  = '%s/ko_C_stats.txt' % op_dir

    ############################################## Read in KEGG DB files ###############################################

    print('Read in KEGG DB files')
    As_description_dict = {}
    Bs_description_dict = {}
    Cs_description_dict = {}
    Ds_description_dict = {}
    D2ABCD_dict = {}
    current_A = ''
    current_B = ''
    current_C = ''
    for each_line in open(db_file):
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

    ABCD_description_dict = {}
    for each_a in As_description_dict:
        ABCD_description_dict[each_a] = As_description_dict[each_a]
    for each_b in Bs_description_dict:
        ABCD_description_dict[each_b] = Bs_description_dict[each_b]
    for each_c in Cs_description_dict:
        ABCD_description_dict[each_c] = Cs_description_dict[each_c]
    for each_d in Ds_description_dict:
        ABCD_description_dict[each_d] = Ds_description_dict[each_d]

    ########################################################################################################################

    input_ko_id_set = set()
    for each_ko in open(ko_txt):
        ko_id = each_ko.strip().split()[0]
        if ko_id.startswith('K'):
            input_ko_id_set.add(ko_id)

    ko_a_stats_dict = dict()
    ko_b_stats_dict = dict()
    ko_c_stats_dict = dict()
    for ko_id in input_ko_id_set:
        ko_abcd_list = D2ABCD_dict.get(ko_id, [])
        for ko_abcd in ko_abcd_list:
            ko_abcd_split = ko_abcd.split('|')
            ko_a = ko_abcd_split[0][2:]
            ko_b = ko_abcd_split[1][2:]
            ko_c = ko_abcd_split[2][2:]
            if ko_a not in ko_a_stats_dict:
                ko_a_stats_dict[ko_a] = set()
            if ko_b not in ko_b_stats_dict:
                ko_b_stats_dict[ko_b] = set()
            if ko_c not in ko_c_stats_dict:
                ko_c_stats_dict[ko_c] = set()
            ko_a_stats_dict[ko_a].add(ko_id)
            ko_b_stats_dict[ko_b].add(ko_id)
            ko_c_stats_dict[ko_c].add(ko_id)

    ko_stats_to_txt(ko_a_stats_dict, ABCD_description_dict, op_a_txt, op_a_stats_txt)
    ko_stats_to_txt(ko_b_stats_dict, ABCD_description_dict, op_b_txt, op_b_stats_txt)
    ko_stats_to_txt(ko_c_stats_dict, ABCD_description_dict, op_c_txt, op_c_stats_txt)

    print('Done!')


if __name__ == "__main__":

    stats_ko_parser = argparse.ArgumentParser(usage=stats_ko_usage)
    stats_ko_parser.add_argument('-ko', required=True,                          help='ko.txt')
    stats_ko_parser.add_argument('-db', required=True,                          help='ko00001.keg')
    stats_ko_parser.add_argument('-o',  required=True,                          help='output directory')
    stats_ko_parser.add_argument('-f',  required=False, action="store_true",    help='force overwrite')
    args = vars(stats_ko_parser.parse_args())
    stats_ko(args)
