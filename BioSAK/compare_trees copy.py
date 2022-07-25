import os
import argparse
from ete3 import Tree
from BioSAK.BioSAK_config import config_dict

compare_trees_usage = '''
============== compare_trees example command ==============

BioSAK compare_trees -t1 tree_1.newick -t2 tree_2.newick

===========================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def compare_trees(args):

    compare_trees_R = config_dict['compare_trees_R']
    tree_file_1 = args['t1']
    tree_file_2 = args['t2']

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
        print('No leaves shared by t1 and t2, program exited!')
        exit()

    elif len(tree1_leaf_list) == len(tree2_leaf_list) == len(shared_leaves):
        compare_trees_cmd = 'Rscript %s -a %s -b %s' % (compare_trees_R, tree_file_1, tree_file_2)
        os.system(compare_trees_cmd)

    elif (len(shared_leaves) != len(tree1_leaf_list)) or (len(shared_leaves) != len(tree2_leaf_list)):
        print('Different leaves were found in t1 (%s) and t2 (%s), will perform mantel test based on shared leaves (%s)' % (
                len(tree1_leaf_list), len(tree2_leaf_list), len(shared_leaves)))

        tree1_path, tree1_basename, tree1_extension = sep_path_basename_ext(tree_file_1)
        tree2_path, tree2_basename, tree2_extension = sep_path_basename_ext(tree_file_2)

        # write out shared leaves
        shared_leaves_txt = '%s_vs_%s_shared_leaves.txt' % (tree1_basename, tree2_basename)
        shared_leaves_txt_handle = open(shared_leaves_txt, 'w')
        for each_shared_leaf in shared_leaves:
            shared_leaves_txt_handle.write(each_shared_leaf + '\n')
        shared_leaves_txt_handle.close()

        # subset_tree
        t1_subset     = '%s_vs_%s_%s_subset%s' % (tree1_basename, tree2_basename, tree1_basename, tree1_extension)
        t2_subset     = '%s_vs_%s_%s_subset%s' % (tree1_basename, tree2_basename, tree2_basename, tree2_extension)
        subset_cmd_t1 = 'BioSAK subset_tree -tree %s -taxon %s -out %s' % (tree_file_1, shared_leaves_txt, t1_subset)
        subset_cmd_t2 = 'BioSAK subset_tree -tree %s -taxon %s -out %s' % (tree_file_2, shared_leaves_txt, t2_subset)
        os.system(subset_cmd_t1)
        os.system(subset_cmd_t2)

        compare_trees_cmd = 'Rscript %s -a %s -b %s' % (compare_trees_R, t1_subset, t2_subset)
        os.system(compare_trees_cmd)


if __name__ == '__main__':

    compare_trees_parser = argparse.ArgumentParser(usage=compare_trees_usage)
    compare_trees_parser.add_argument('-t1', required=True, help='tree 1')
    compare_trees_parser.add_argument('-t2', required=True, help='tree 2')
    args = vars(compare_trees_parser.parse_args())

    compare_trees(args)
