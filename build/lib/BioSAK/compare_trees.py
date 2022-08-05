import os
import glob
import argparse
from ete3 import Tree
import multiprocessing as mp
from BioSAK.BioSAK_config import config_dict

compare_trees_usage = '''
================ compare_trees example command ================

module load R/3.5.3
BioSAK compare_trees -t1 tree_1.newick -t2 tree_2.newick
BioSAK compare_trees -t1 tree_dir -t2 tree_dir -tx newick -dm -t 12

===============================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def check_numeric(str_in):
    is_numeric = True
    try:
        x = float(str_in)
    except ValueError:
        is_numeric = False

    return is_numeric


def parse_mantel_stats(mantel_stats_txt):

    mantel_similarity = 'na'
    for each_line in open(mantel_stats_txt):
        if 'Mantel statistic r: ' in each_line:
            mantel_similarity = each_line.strip().split('Mantel statistic r: ')[1]
    return mantel_similarity


def get_matrix(query_tree_list, subject_tree_list, mantel_stats_dir, write_out_dm, output_matrix, output_matrix_distance):

    header_line_str = '\t' + '\t'.join(subject_tree_list) + '\n'

    output_matrix_handle = open(output_matrix, 'w')
    output_matrix_handle.write(header_line_str)
    distance_lol = []
    for each_qt in query_tree_list:

        current_qt_mantel_stats_value_list = [each_qt]
        for each_st in subject_tree_list:

            qt_vs_st_mantel_stats = '%s/%s_vs_%s_mantel_stats.txt' % (mantel_stats_dir, each_qt, each_st)
            st_vs_qt_mantel_stats = '%s/%s_vs_%s_mantel_stats.txt' % (mantel_stats_dir, each_st, each_qt)

            tree_similarity = 'na'
            if os.path.isfile(qt_vs_st_mantel_stats) is True:
                tree_similarity = parse_mantel_stats(qt_vs_st_mantel_stats)
            if os.path.isfile(st_vs_qt_mantel_stats) is True:
                tree_similarity = parse_mantel_stats(st_vs_qt_mantel_stats)

            current_qt_mantel_stats_value_list.append(tree_similarity)

        current_qt_mantel_stats_value_list_distance = [each_qt]
        for each_value in current_qt_mantel_stats_value_list[1:]:
            if check_numeric(each_value) is True:
                in_distance = 1 - float(each_value)
                in_distance = float("{0:.4f}".format(in_distance))
                if in_distance == 0:
                    in_distance = '0'
                current_qt_mantel_stats_value_list_distance.append(str(in_distance))
            else:
                current_qt_mantel_stats_value_list_distance.append('na')

        distance_lol.append(current_qt_mantel_stats_value_list_distance)
        current_qt_mantel_stats_value_str = '\t'.join(current_qt_mantel_stats_value_list)
        output_matrix_handle.write(current_qt_mantel_stats_value_str + '\n')
    output_matrix_handle.close()

    # write out distance matrix
    if write_out_dm is True:
        output_matrix_distance_handle = open(output_matrix_distance, 'w')
        output_matrix_distance_handle.write(header_line_str)
        for each_list in distance_lol:
            output_matrix_distance_handle.write('\t'.join(each_list) + '\n')
        output_matrix_distance_handle.close()


def compare_trees_worker(arg_list):

    compare_trees_R = arg_list[0]
    tree_file_1     = arg_list[1]
    tree_file_2     = arg_list[2]
    keep_tmp_file   = arg_list[3]

    tree1_path, tree1_basename, tree1_extension = sep_path_basename_ext(tree_file_1)
    tree2_path, tree2_basename, tree2_extension = sep_path_basename_ext(tree_file_2)

    op_stats = '%s_vs_%s_mantel_stats.txt' % (tree1_basename, tree2_basename)

    t1 = Tree(tree_file_1, format=1)
    t2 = Tree(tree_file_2, format=1)

    tree1_leaf_list = []
    for leaf1 in t1:
        tree1_leaf_list.append(leaf1.name)

    tree2_leaf_list = []
    for leaf2 in t2:
        tree2_leaf_list.append(leaf2.name)

    shared_leaves = set(tree1_leaf_list).intersection(tree2_leaf_list)
    if len(shared_leaves) == 0:
        print('No leaves shared between %s and %s, calculation skipped!' % (tree1_basename, tree2_basename))
        #exit()

    elif len(tree1_leaf_list) == len(tree2_leaf_list) == len(shared_leaves):
        compare_trees_cmd = 'Rscript %s -a %s -b %s > %s' % (compare_trees_R, tree_file_1, tree_file_2, op_stats)
        os.system(compare_trees_cmd)

    elif (len(shared_leaves) != len(tree1_leaf_list)) or (len(shared_leaves) != len(tree2_leaf_list)):
        print('Performing Mantel test based on %s leaves shared by %s (%s) and %s (%s)' % (len(shared_leaves), tree1_basename, len(tree1_leaf_list), tree2_basename, len(tree2_leaf_list)))

        # write out shared leaves
        shared_leaves_txt = '%s_vs_%s_shared_leaves.txt' % (tree1_basename, tree2_basename)
        shared_leaves_txt_handle = open(shared_leaves_txt, 'w')
        for each_shared_leaf in shared_leaves:
            shared_leaves_txt_handle.write(each_shared_leaf + '\n')
        shared_leaves_txt_handle.close()

        # subset_tree
        t1_subset     = '%s_vs_%s_%s_subset%s' % (tree1_basename, tree2_basename, tree1_basename, tree1_extension)
        t2_subset     = '%s_vs_%s_%s_subset%s' % (tree1_basename, tree2_basename, tree2_basename, tree2_extension)
        subset_cmd_t1 = 'BioSAK subset_tree -tree %s -taxon %s -out %s -q' % (tree_file_1, shared_leaves_txt, t1_subset)
        subset_cmd_t2 = 'BioSAK subset_tree -tree %s -taxon %s -out %s -q' % (tree_file_2, shared_leaves_txt, t2_subset)
        os.system(subset_cmd_t1)
        os.system(subset_cmd_t2)

        compare_trees_cmd = 'Rscript %s -a %s -b %s > %s' % (compare_trees_R, t1_subset, t2_subset, op_stats)
        os.system(compare_trees_cmd)

        if keep_tmp_file is False:
            os.system('rm %s' % shared_leaves_txt)
            os.system('rm %s' % t1_subset)
            os.system('rm %s' % t2_subset)


def compare_trees(args):

    compare_trees_R = config_dict['compare_trees_R']
    #output_prefix   = args['p']
    tree_file_1     = args['t1']
    tree_file_2     = args['t2']
    tree_file_ext   = args['tx']
    export_dm       = args['dm']
    num_threads     = args['t']
    keep_tmp        = args['tmp']

    query_tree_list = []
    if os.path.isfile(tree_file_1):
        query_tree_list = [tree_file_1]
    elif os.path.isdir(tree_file_1):
        query_tree_re = '%s/*.%s' % (tree_file_1, tree_file_ext)
        query_tree_list = glob.glob(query_tree_re)

    subject_tree_list = []
    if os.path.isfile(tree_file_2):
        subject_tree_list = [tree_file_2]
    elif os.path.isdir(tree_file_2):
        subject_tree_re = '%s/*.%s' % (tree_file_2, tree_file_ext)
        subject_tree_list = glob.glob(subject_tree_re)

    # prepare arg list for compare_trees_worker
    to_be_calculated_set = set()
    list_for_compare_trees_worker = []
    for each_query_tree in query_tree_list:
        for each_subject_tree in subject_tree_list:

            tree_1_vs_2 = '%s_vs_%s' % (each_query_tree, each_subject_tree)
            tree_2_vs_1 = '%s_vs_%s' % (each_subject_tree, each_query_tree)

            if tree_1_vs_2 not in to_be_calculated_set:
                list_for_compare_trees_worker.append([compare_trees_R, each_query_tree, each_subject_tree, keep_tmp])
                to_be_calculated_set.add(tree_1_vs_2)
                to_be_calculated_set.add(tree_2_vs_1)

    print('Total pairs of trees to compare: %s' % len(list_for_compare_trees_worker))

    # compare trees with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(compare_trees_worker, list_for_compare_trees_worker)
    pool.close()
    pool.join()

    # get matrix
    output_matrix_similarity    = 'Matrix_similarity.txt'
    output_matrix_distance      = 'Matrix_distance.txt'
    query_tree_list_basename = []
    for each_q_tree in query_tree_list:
        q_tree_path, q_tree_basename, q_tree_ext = sep_path_basename_ext(each_q_tree)
        query_tree_list_basename.append(q_tree_basename)

    subject_tree_list_basename = []
    for each_s_tree in subject_tree_list:
        s_tree_path, s_tree_basename, s_tree_ext = sep_path_basename_ext(each_s_tree)
        subject_tree_list_basename.append(s_tree_basename)

    get_matrix(sorted(query_tree_list_basename), sorted(subject_tree_list_basename), '.', export_dm, output_matrix_similarity, output_matrix_distance)

    # final report
    if export_dm is True:
        print('Data matrix exported to: %s and %s' % (output_matrix_similarity, output_matrix_distance))
    else:
        print('Data matrix exported to: %s' % output_matrix_similarity)

    print('Done!')


if __name__ == '__main__':

    compare_trees_parser = argparse.ArgumentParser(usage=compare_trees_usage)
    #compare_trees_parser.add_argument('-p',   required=True,                       help='output prefix')
    compare_trees_parser.add_argument('-t1',  required=True,                       help='tree (folder) 1')
    compare_trees_parser.add_argument('-t2',  required=True,                       help='tree (folder) 2')
    compare_trees_parser.add_argument('-tx',  required=False, default='newick',    help='extention of tree files, default: newick')
    compare_trees_parser.add_argument('-dm',  required=False, action="store_true", help='export distance-alike matrix, obtained by subtract the similarity value from 1')
    compare_trees_parser.add_argument('-t',   required=False, type=int, default=1, help='number of threads')
    compare_trees_parser.add_argument('-tmp', required=False, action="store_true", help='keep tmp files')
    args = vars(compare_trees_parser.parse_args())
    compare_trees(args)
