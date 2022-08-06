import argparse
from ete3 import Tree


rename_leaves_usage = '''
==================== rename_leaves example commands ====================

BioSAK rename_leaves -tree input.tree -txt rename.txt -out output.tree

# format of rename.txt (tab separated)
leaf_1  leaf_1_new_name
leaf_2  leaf_2_new_name

========================================================================
'''


def rename_leaves(args):

    tree_file_in      = args['tree']
    rename_file       = args['txt']
    tree_format       = args['fmt']
    tree_file_out     = args['out']

    mag_rename_dict = {}
    for each_mag in open(rename_file):
        each_mag_split = each_mag.strip().split('\t')
        before_rename = each_mag_split[0]
        after_rename = each_mag_split[1]
        mag_rename_dict[after_rename] = before_rename

    t = Tree(tree_file_in, format=tree_format)

    input_tree_leaf_name_list = []
    for leaf in t:
        input_tree_leaf_name_list.append(leaf.name)

    leaves_with_new_name = 0
    leaves_without_new_name = 0
    for each_raw_name in input_tree_leaf_name_list:
        if each_raw_name in mag_rename_dict:
            leaves_with_new_name += 1
        else:
            leaves_without_new_name += 1

    if leaves_with_new_name == 0:
        print('No leaf on input tree found in rename file, please double check!')
        exit()
    elif leaves_without_new_name > 0:
        print('Found %s leaves in input tree, %s of them were found in the rename file' % (len(input_tree_leaf_name_list), leaves_with_new_name))

    for leaf in t:
        leaf_name_new = mag_rename_dict.get(leaf.name, leaf.name)
        leaf.name = leaf_name_new
    t.write(format=tree_format, outfile=tree_file_out)


if __name__ == '__main__':

    rename_leaves_parser = argparse.ArgumentParser()
    rename_leaves_parser.add_argument('-tree', required=True,             help='input tree')
    rename_leaves_parser.add_argument('-txt',  required=True,             help='rename file')
    rename_leaves_parser.add_argument('-fmt',  required=False, default=1, help='tree format, default: 1')
    rename_leaves_parser.add_argument('-out',  required=True,             help='output tree')
    args = vars(rename_leaves_parser.parse_args())
    rename_leaves(args)
