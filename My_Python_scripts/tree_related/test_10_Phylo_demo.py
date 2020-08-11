import copy
from Bio import Phylo


def get_child_nodes(tree, node_id):
    for clade in tree.find_clades():
        if str(clade.name) == node_id:
            all_children_nodes = clade.get_terminals() + clade.get_nonterminals()

            all_children_nodes_self_excluded = []
            for nonleaf_node in all_children_nodes:
                if str(nonleaf_node.name) != node_id:
                    all_children_nodes_self_excluded.append(str(nonleaf_node.name))

            return all_children_nodes_self_excluded


# not perfect
def subset_tree(tree, taxon_to_keep):

    # add tmp value to None clade
    name_index = 0
    for clade in tree.find_clades():
        if clade.name == None:
            clade.name = 'TmpAdd_%s' % name_index
        name_index += 1

    for clade in tree.find_clades():
        if str(clade.name) in taxon_to_keep:
            clade_child_nodes = get_child_nodes(tree, str(clade.name))
            for child_node in clade_child_nodes:
                tree.collapse(child_node)

    all_leaf_clades = tree.get_terminals()
    all_leaf_clades_id = [str(i.name) for i in all_leaf_clades]
    all_leaf_clades_id_to_collapse = []
    for each_leaf in all_leaf_clades_id:
        if each_leaf not in taxon_to_keep:
            all_leaf_clades_id_to_collapse.append(each_leaf)

    for leaf_to_collapse in all_leaf_clades_id_to_collapse:
        tree.collapse(leaf_to_collapse)

    return tree


######################################################### demo #########################################################

tree_file_in = '/Users/songweizhi/Desktop/tree/example.tree'
tree_in = Phylo.read(tree_file_in, 'newick')
tree_new = copy.deepcopy(tree_in)


#tree_new.get_terminals('X', order='postorder')
# , order='preorder' postorder level
#name_list = ['A']
#for name in name_list:
    #tree_new.prune(name)
    #tree_new.collapse(name)


#print(tree_new.collapse_all('X'))
#tree_new.collapse('A')
#tree_new.collapse('B')
#tree_new.collapse('I')
#tree_new.prune('I')

#print(tree_new.find_any('X'))





# for clade in tree_new.find_clades():
#     if str(clade.name) in taxon_to_keep:
#         clade_child_nodes = get_child_nodes(tree_new, str(clade.name))
#
#         print(clade_child_nodes)
#         # for child_node in clade_child_nodes:
#         #
#         #     print(child_node)
#         #     #tree.collapse(child_node)

#
# taxon_to_keep = ['H', 'I']
# tree_new = subset_tree(tree_new, taxon_to_keep)
# tree_new = subset_tree(tree_new, taxon_to_keep)
# tree_new = subset_tree(tree_new, taxon_to_keep)
# tree_new = subset_tree(tree_new, taxon_to_keep)
# tree_new = subset_tree(tree_new, taxon_to_keep)
#

# to compare
print('\n########## For comparison ##########\n')
print(tree_in)
print('\n')
print(tree_new)
print('\n')
Phylo.draw_ascii(tree_in)
print('\n')
Phylo.draw_ascii(tree_new)
