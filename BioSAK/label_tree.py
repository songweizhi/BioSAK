import os
import argparse
from BioSAK.BioSAK_config import config_dict


label_tree_usage = '''
======================== label_tree example commands ========================

module load R

# label tree with customized grouping file
BioSAK label_tree -tree NorthSea.tree -label labels.txt 

# label tree by taxonomic classification at phylum and class levels
BioSAK label_tree -tree NorthSea.tree -taxon GTDB_output.tsv -rank p
BioSAK label_tree -tree NorthSea.tree -taxon GTDB_output.tsv -rank c

# label file format:
label_A,tree_leaf_1
label_B,tree_leaf_2
label_B,tree_leaf_3
label_C,tree_leaf_4

=============================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def label_tree(args, config_dict):

    tree_in      = args['tree']
    label_file   = args['label']
    leaf_taxon   = args['taxon']
    taxon_rank   = args['rank']
    label_tree_R = config_dict['label_tree_R']

    if (label_file is not None) and (leaf_taxon is None) and (taxon_rank is None):
        label_tree_cmd = 'Rscript %s -t %s -g %s' % (label_tree_R, tree_in, label_file)
        os.system(label_tree_cmd)

    elif (label_file is None) and (leaf_taxon is not None) and (taxon_rank is not None):

        # define tmp file name
        tree_file_path, tree_file_basename, tree_file_extension = sep_path_basename_ext(tree_in)
        taxon_grouping = '%s/%s_%s.txt' % (tree_file_path, tree_file_basename, taxon_rank)

        # read GTDB output into dict
        taxon_assignment_dict = {}
        for each_genome in open(leaf_taxon):
            if not each_genome.startswith('user_genome'):
                each_split = each_genome.strip().split('\t')
                bin_name = each_split[0]

                assignment_full = []
                if len(each_split) == 1:
                    assignment_full = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
                elif (len(each_split) > 1) and (';' in each_split[1]):
                    assignment = each_split[1].split(';')
                    if len(assignment) == 7:
                        assignment_full = assignment
                    if len(assignment) == 6:
                        assignment_full = assignment + ['s__']
                    if len(assignment) == 5:
                        assignment_full = assignment + ['g__', 's__']
                    if len(assignment) == 4:
                        assignment_full = assignment + ['f__', 'g__', 's__']
                    if len(assignment) == 3:
                        assignment_full = assignment + ['o__', 'f__', 'g__', 's__']
                    if len(assignment) == 2:
                        assignment_full = assignment + ['c__', 'o__', 'f__', 'g__', 's__']

                elif (len(each_split) > 1) and (';' not in each_split[1]):
                    assignment_full = [each_split[1]] + ['p__', 'c__', 'o__', 'f__', 'g__', 's__']

                # store in dict
                taxon_assignment_dict[bin_name] = assignment_full

        # get all identified taxon at defined ranks
        rank_to_position_dict = {'d': 0, 'p': 1, 'c': 2, 'o': 3, 'f': 4, 'g': 5, 's': 6}
        specified_rank_pos = rank_to_position_dict[taxon_rank]

        taxon_grouping_handle = open(taxon_grouping, 'w')
        for each_TaxonAssign in taxon_assignment_dict:
            specified_rank_id = taxon_assignment_dict[each_TaxonAssign][specified_rank_pos]
            taxon_grouping_handle.write('%s,%s\n' % (specified_rank_id, each_TaxonAssign))
        taxon_grouping_handle.close()

        # run R script
        label_tree_cmd = 'Rscript %s -t %s -g %s' % (label_tree_R, tree_in, taxon_grouping)
        os.system(label_tree_cmd)

    else:
        print('Please provide either a customized label file or the taxonomy info of tree leaves together with a taxonomic rank')
        print('Program exited!')
        exit()


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser(usage=label_tree_usage)

    parser.add_argument('-tree',      required=True,                 help='tree file in newick format')
    parser.add_argument('-label',     required=False,  default=None, help='label file (label,leaf)')
    parser.add_argument('-taxon',     required=False,  default=None, help='taxonomic classification')
    parser.add_argument('-rank',      required=False,  default=None, help='taxonomic rank to label')

    args = vars(parser.parse_args())
    label_tree(args, config_dict)
