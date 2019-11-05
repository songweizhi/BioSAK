
newick_file = '/Users/songweizhi/Desktop/human_gut_species_tree.newick'
species_list = ['AKV', 'BGC', 'BAD']

subtree_file = '/Users/songweizhi/Desktop/human_gut_species_tree_subset.newick'


###################################################### Biopython #######################################################

from Bio import Phylo

# tree = Phylo.read(newick_file, 'newick')
# common_ancestor = tree.common_ancestor(species_list)
# Phylo.draw_ascii(common_ancestor)


######################################################### ETE3 #########################################################

from ete3 import Tree


def subset_tree(tree_file_in, leaf_node_list, tree_file_out):
    tree_in = Tree(tree_file_in, format=0)
    tree_in.prune(leaf_node_list, preserve_branch_length=True)
    tree_in.write(format=0, outfile=tree_file_out)


subset_tree(newick_file, species_list, subtree_file)

