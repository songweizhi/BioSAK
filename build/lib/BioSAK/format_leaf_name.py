import os
import argparse
from ete3 import Tree


format_leaf_name_usage = '''
========= FLN (Format Leaf Name) example commands =========

BioSAK FLN -i input.tree -o output.tree -s2u -nsqm -ndqm
BioSAK FLN -i input.tree -o output.tree -ns -nsqm

===========================================================
'''


def format_leaf_name(args):

    tree_file_in             = args['i']
    tree_format              = args['fmt']
    tree_file_out            = args['o']
    no_space                 = args['ns']
    space_to_underscore      = args['s2u']
    no_single_quotation_mark = args['nsqm']
    no_double_quotation_mark = args['ndqm']

    if os.path.isfile(tree_file_in) is False:
        print('Tree file not found, program exited!')
        exit()

    if (no_space is True) and (space_to_underscore is True):
        print('Two actions (-ns and -s2u) specified to spaces in tree leaves, program exited!')
        exit()

    t = Tree(tree_file_in, format=tree_format)

    # get rename dict
    mag_rename_dict = dict()
    for leaf in t:
        leaf_name = leaf.name
        if space_to_underscore is True:
            leaf_name = leaf_name.replace(' ', '_')
        if no_space is True:
            leaf_name = leaf_name.replace(' ', '')
        if no_single_quotation_mark is True:
            leaf_name = leaf_name.replace("'", '')
        if no_double_quotation_mark is True:
            leaf_name = leaf_name.replace('"', '')
        mag_rename_dict[leaf.name] = leaf_name

    for leaf in t:
        leaf_name = leaf.name
        leaf_name_new = mag_rename_dict[leaf_name]
        leaf.name = leaf_name_new
    t.write(format=tree_format, outfile=tree_file_out)

    print('Done!')


if __name__ == '__main__':

    format_leaf_name_parser = argparse.ArgumentParser()
    format_leaf_name_parser.add_argument('-i',                  required=True,                          help='input tree')
    format_leaf_name_parser.add_argument('-fmt',                required=False, default=1,              help='tree format, default: 1')
    format_leaf_name_parser.add_argument('-o',                  required=True,                          help='output tree')
    format_leaf_name_parser.add_argument('-s2u',                required=False, action="store_true",    help='change space in tree leaves to underscore')
    format_leaf_name_parser.add_argument('-ns',                 required=False, action="store_true",    help='remove space from leaf names')
    format_leaf_name_parser.add_argument('-nsqm',               required=False, action="store_true",    help='remove single quotation marks from leaf names')
    format_leaf_name_parser.add_argument('-ndqm',               required=False, action="store_true",    help='remove double quotation marks from leaf names')
    args = vars(format_leaf_name_parser.parse_args())
    format_leaf_name(args)
