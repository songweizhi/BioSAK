import os
import argparse


cdb_usage = '''
===== cdb example commands =====

BioSAK cdb -i Cdb.csv

================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def cdb(args):

    Cdb_file = args['i']
    cdb_name, cdb_path, cdb_base, cdb_ext = sep_path_basename_ext(Cdb_file)

    op_file_1 = '%s/%s_reformatted_1.%s' % (cdb_path, cdb_base, cdb_ext)
    op_file_2 = '%s/%s_reformatted_2.%s' % (cdb_path, cdb_base, cdb_ext)

    bin_to_cluster_dict = {}
    cluster_to_bin_dict = {}
    obtained_clusters = set()
    col_index = dict()
    line_num_index = 0
    for each_line in open(Cdb_file):
        line_num_index += 1
        line_split = each_line.strip().split(',')
        if line_num_index == 1:
            col_index = {key: i for i, key in enumerate(line_split)}
        else:
            bin_id = line_split[col_index['genome']]
            secondary_cluster = line_split[col_index['secondary_cluster']]
            obtained_clusters.add(secondary_cluster)
            bin_to_cluster_dict[bin_id] = secondary_cluster
            if secondary_cluster not in cluster_to_bin_dict:
                cluster_to_bin_dict[secondary_cluster] = [bin_id]
            else:
                cluster_to_bin_dict[secondary_cluster].append(bin_id)

    obtained_clusters_list = sorted([i for i in obtained_clusters])

    # write out 1
    op_file_1_handle = open(op_file_1, 'w')
    for j in obtained_clusters_list:
        op_file_1_handle.write('cluster_%s\t%s\n' % (j, '\t'.join(cluster_to_bin_dict[j])))
    op_file_1_handle.close()

    # write out 2
    op_file_2_handle = open(op_file_2, 'w')
    for j in obtained_clusters_list:
        cluster_gnm_set = cluster_to_bin_dict[j]
        for each_g in cluster_gnm_set:
            op_file_2_handle.write('%s\t%s\n' % (each_g, j))
    op_file_2_handle.close()


if __name__ == '__main__':
    cdb_parser = argparse.ArgumentParser(usage=cdb_usage)
    cdb_parser.add_argument('-i',  required=True, help='input file, Cdb.csv')
    args = vars(cdb_parser.parse_args())
    cdb(args)

