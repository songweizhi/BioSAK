from ete3 import NCBITaxa, Tree, TreeStyle, NodeStyle, TextFace


def plot_tree(tree_newick, tree_title):

    tree = Tree(tree_newick, format = 1)
    # set tree parameters
    ts = TreeStyle()
    ts.mode = "r"  # tree model: 'r' for rectangular, 'c' for circular
    ts.show_leaf_name = 0
    # set tree title text parameters
    ts.title.add_face(TextFace(tree_title,
                               fsize = 8,
                               fgcolor = 'black',
                               ftype = 'Arial',
                               tight_text = False),
                      column = 0)  # tree title text setting
    # set layout parameters
    ts.rotation = 0  # from 0 to 360
    ts.show_scale = False
    ts.margin_top = 10  # top tree image margin
    ts.margin_bottom = 10  # bottom tree image margin
    ts.margin_left = 10  # left tree image margin
    ts.margin_right = 10  # right tree image margin
    ts.show_border = False  # set tree image border
    ts.branch_vertical_margin = 3  # 3 pixels between adjancent branches

    # set tree node style
    for each_node in tree.traverse():
        # leaf node parameters
        if each_node.is_leaf():
            ns = NodeStyle()
            ns["shape"] = "circle"  # dot shape: circle, square or sphere
            ns["size"] = 0  # dot size
            ns['hz_line_width'] = 0.5  # branch line width
            ns['vt_line_width'] = 0.5  # branch line width
            ns['hz_line_type'] = 0  # branch line type: 0 for solid, 1 for dashed, 2 for dotted
            ns['vt_line_type'] = 0  # branch line type
            ns["fgcolor"] = "blue"  # the dot setting
            each_node.add_face(TextFace(each_node.name,
                                        fsize = 5,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 0,
                               position = 'branch-right')  # leaf node the node name text setting

            each_node.set_style(ns)

        # non-leaf node parameters
        else:
            nlns = NodeStyle()
            nlns["size"] = 0  # dot size
            #nlns["rotation"] = 45
            each_node.add_face(TextFace(each_node.name[:12],

                                        fsize = 4,
                                        fgcolor = 'black',
                                        tight_text = False,
                                        bold = False),
                               column = 5,
                               position = 'branch-top')  # non-leaf node name text setting)

            each_node.set_style(nlns)

    tree.render('/Users/songweizhi/Desktop/Tree_' + tree_title + '.png', w = 900, units = "px", tree_style = ts)  # set figures size


taxon_id_assignment_file = '/Users/songweizhi/Desktop/taxon_assignment_168.txt'
taxdump_file = '/Users/songweizhi/Desktop/taxdump_customized.tar.gz'

# get taxon_id_list
taxon_ids = open(taxon_id_assignment_file)
taxon_id_list = []
for each in taxon_ids:
    each_split = each.strip().split('\t')
    taxon_id = each_split[1]
    taxon_id_list.append(taxon_id)

# print(taxon_id_list)
# print(len(taxon_id_list))

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database(taxdump_file)

tree_phylo = ncbi.get_topology(taxon_id_list, intermediate_nodes=1)

tree_phylo_newick = tree_phylo.write(format=8)

for each_node in tree_phylo.traverse():
    node_name_list = []
    node_name_list.append(each_node.name)
    if node_name_list == ['']:
        pass
    else:
        if each_node.is_leaf():
            # change bin id to bin name
            each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]
            # add group information to bin name

        # for non-leaf node name, only display the first 12 characters
        else:
            if len(ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]) <= 12:
                each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)]
            else:
                each_node.name = ncbi.get_taxid_translator(node_name_list)[int(each_node.name)][0:12] + '.'

print(tree_phylo)
tree_phylo_newick_2 = tree_phylo.write(format=8)

print(tree_phylo_newick_2)

plot_tree(tree_phylo_newick, '168 bins_id')
plot_tree(tree_phylo_newick_2, '168 bins_name')







