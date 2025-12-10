import math
import os.path
import random
import argparse
import seaborn as sns


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']
    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']
    else:
        color_num_each = math.ceil(color_num/8) + 2
        color_list_1 = sns.color_palette('Blues', n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',  n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()
        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    random.shuffle(color_list_to_return_sorted)

    return color_list_to_return_sorted


def cdb_to_itol_piechart(dereplicated_gnm_id_txt, cdb_tsv, metadata_txt, gnm_id_col, interested_grp_col, group_color_txt, op_txt, op_txt_itol):

    grp_dict = dict()
    col_index = dict()
    for each_gnm in open(metadata_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        if each_gnm.startswith('%s\t' % gnm_id_col):
            col_index = {key: i for i, key in enumerate(each_gnm_split)}
        else:
            gnm_id           = each_gnm_split[col_index[gnm_id_col]]
            interested_grp   = each_gnm_split[col_index[interested_grp_col]]
            grp_dict[gnm_id] = interested_grp

    dereplicated_gnm_set = set()
    for each_gnm in open(dereplicated_gnm_id_txt):
        dereplicated_gnm_set.add(each_gnm.strip().split()[0])

    bin_to_cluster_dict = {}
    cluster_to_bin_dict = {}
    obtained_clusters = set()
    col_index = dict()
    line_num_index = 0
    for each_line in open(cdb_tsv):
        line_num_index += 1
        line_split = each_line.strip().split(',')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            gnm_name = line_split[col_index['genome']]
            gnm_id   = '.'.join(gnm_name.split('.')[:-1])
            secondary_cluster = line_split[col_index['secondary_cluster']]
            obtained_clusters.add(secondary_cluster)
            bin_to_cluster_dict[gnm_id] = secondary_cluster
            if secondary_cluster not in cluster_to_bin_dict:
                cluster_to_bin_dict[secondary_cluster] = [gnm_id]
            else:
                cluster_to_bin_dict[secondary_cluster].append(gnm_id)

    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('Cluster\tRepresentative\tMembers\n')
    grp_dod = dict()
    all_grp_set = set()
    for each_gnm in dereplicated_gnm_set:
        c_id = bin_to_cluster_dict[each_gnm]
        c_gnm_list_sorted = sorted(list(cluster_to_bin_dict[c_id]))
        grp_list = [grp_dict[i] for i in c_gnm_list_sorted]
        current_grp_dict = dict()
        for grp in grp_list:
            if grp not in current_grp_dict:
                current_grp_dict[grp] = 0
            current_grp_dict[grp] += 1
            all_grp_set.add(grp)
        grp_dod[each_gnm] = current_grp_dict
        op_txt_handle.write('%s\t%s\t%s\n' % (c_id, each_gnm, ','.join(c_gnm_list_sorted)))
    op_txt_handle.close()

    all_grp_list_sorted = sorted(list(all_grp_set))
    grp_color_dict = dict()
    if group_color_txt is not None:
        for each_g in open(group_color_txt):
            each_g_split = each_g.strip().split('\t')
            grp_color_dict[each_g_split[0]] = each_g_split[1]
    else:
        grp_color_list = get_color_list(len(all_grp_list_sorted))
        for (grp, color) in zip(all_grp_list_sorted, grp_color_list):
            grp_color_dict[grp] = color
    grp_color_list = [grp_color_dict.get(i, '#000000') for i in all_grp_list_sorted]

    # write out itol file
    op_txt_itol_handle = open(op_txt_itol, 'w')
    op_txt_itol_handle.write('DATASET_PIECHART\n')
    op_txt_itol_handle.write('SEPARATOR TAB\n')
    op_txt_itol_handle.write('DATASET_LABEL\tPieChart\n')
    op_txt_itol_handle.write('FIELD_COLORS\t%s\n' % '\t'.join([grp_color_dict.get(i, '#000000') for i in all_grp_list_sorted]))
    op_txt_itol_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(all_grp_list_sorted))
    op_txt_itol_handle.write('HEIGHT_FACTOR\t1.8\n')
    op_txt_itol_handle.write('BORDER_WIDTH\t0\n')
    op_txt_itol_handle.write('LEGEND_TITLE\tPieChart\n')
    op_txt_itol_handle.write('LEGEND_SCALE\t1\n')
    op_txt_itol_handle.write('LEGEND_SHAPES\t%s\n' % ('\t'.join(['1']*len(all_grp_list_sorted))))
    op_txt_itol_handle.write('LEGEND_COLORS\t%s\n' % ('\t'.join(grp_color_list)))
    op_txt_itol_handle.write('LEGEND_LABELS\t%s\n' % ('\t'.join(all_grp_list_sorted)))
    op_txt_itol_handle.write('LEGEND_SHAPE_SCALES\t%s\n' % ('\t'.join(['1']*len(all_grp_list_sorted))))
    op_txt_itol_handle.write('\n')
    op_txt_itol_handle.write('\n')
    op_txt_itol_handle.write('DATA\n')
    op_txt_itol_handle.write('#ID\tPosition\tRadius\t%s\n' % ('\t'.join(all_grp_list_sorted)))

    for each_grp in sorted(list(grp_dod.keys())):
        grp_dict = grp_dod[each_grp]
        value_list = []
        for each_g in all_grp_list_sorted:
            g_count = grp_dict.get(each_g, 0)
            value_list.append(g_count)
        value_list_as_str = [str(i) for i in value_list]
        op_txt_itol_handle.write('%s\t-1\t10\t%s\n' % (each_grp, '\t'.join(value_list_as_str)))
    op_txt_itol_handle.close()


dereplicated_gnm_id_txt = '/Users/songweizhi/Desktop/Sponge_r226/03_AOA_genomes_dRep/Cdb_ANI95_dereplicated_genomes.txt'
cdb_tsv                 = '/Users/songweizhi/Desktop/Sponge_r226/03_AOA_genomes_dRep/Cdb_ANI95.csv'
metadata_txt            = '/Users/songweizhi/Desktop/Sponge_r226/AOA_genomes_combined_1221_metadata.txt'
gnm_id_col              = 'Genome'
interested_grp_col      = 'Host_Category'
group_color_txt         = '/Users/songweizhi/Desktop/Sponge_r226/00_Sponge_r220/0_metadata/color_code_habitat.txt'
op_txt                  = '/Users/songweizhi/Desktop/Sponge_r226/03_AOA_genomes_dRep/aaa.csv'
op_txt_itol             = '/Users/songweizhi/Desktop/Sponge_r226/03_AOA_genomes_dRep/aaa_iTOL.csv'


cdb_to_itol_piechart(dereplicated_gnm_id_txt, cdb_tsv, metadata_txt, gnm_id_col, interested_grp_col, group_color_txt, op_txt, op_txt_itol)
